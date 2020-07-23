/****************************************************************************
 *
 *  CalibPPS/Alignment.plugins/PPSAlignmentWorker.cc
 *
 *  Description : PPS Alignment DQM worker
 *
 * Authors:
 *  - Jan Kašpar
 *  - Mateusz Kocot
 *
 ****************************************************************************/

#include "DQMServices/Core/interface/DQMEDAnalyzer.h"
#include "DQMServices/Core/interface/DQMStore.h"
#include "DQMServices/Core/interface/MonitorElement.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Utilities/interface/ESGetToken.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/CTPPSDetId/interface/TotemRPDetId.h"
#include "DataFormats/CTPPSReco/interface/CTPPSLocalTrackLite.h"
#include "DataFormats/CTPPSReco/interface/CTPPSLocalTrackLiteFwd.h"
#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/CTPPSReco/interface/TotemRPLocalTrack.h"

#include "CondFormats/PPSObjects/interface/PPSAlignmentConfig.h"
#include "CondFormats/DataRecord/interface/PPSAlignmentConfigRcd.h"

#include <vector>
#include <string>
#include <cmath>

#include "TVectorD.h"
#include "TH2D.h"

//----------------------------------------------------------------------------------------------------

class PPSAlignmentWorker : public DQMEDAnalyzer
{
public:
    PPSAlignmentWorker(const edm::ParameterSet &);
    ~PPSAlignmentWorker() override {};

private:
    void bookHistograms(DQMStore::IBooker &, edm::Run const &, edm::EventSetup const &) override;

    void analyze(const edm::Event &, const edm::EventSetup &) override;

    // ------------ structures ------------
    struct Stat
    {
        unsigned int dim;
        double S1;
        std::vector<double> Sv, Svv, Svvv, Svvvv;
        std::vector< std::vector<double> > Sxy, Sxxy, Sxyy, Sxxyy;
        std::vector<std::string> labels;

        Stat();
        Stat(unsigned int dim_);

        void init(unsigned int dim_ = 1);
        void setLabels(const std::vector<std::string> &labels_);

        template<class T>
        void fill(const T &v);
        void fill(double v1, double v2 = 0., double v3 = 0., double v4 = 0., double v5 = 0.);

        std::string qLabel(unsigned int i) const;

        // --- 1D getters ---
        double entries() const;
        double mean(unsigned int i) const;
        double stdDev(unsigned int i) const;
        double meanUnc(unsigned int i) const;
        double stdDevUnc(unsigned int i) const;
        // approximation of GetStdDevUnc valid for Gaussian distributions
        double stdDevUncGauss(unsigned int i) const;

        // --- 2D getters ---
        double covariance(unsigned int i, unsigned int j) const;
        double correlation(unsigned int i, unsigned int j) const;
        double covarianceUnc(unsigned int i, unsigned int j) const;
        double correlationUnc(unsigned int i, unsigned int j) const;
        TMatrixDSym covarianceMatrix() const;

        // --- print methods ---
        void printStat() const;
        void printMeanAndStdDev() const;
        void printCovariance() const;
        void printCorrelation() const;
    };

    struct Profile
    {
        MonitorElement *h = nullptr;
        std::vector<Stat> st;

        MonitorElement *h_entries = nullptr;
        MonitorElement *h_mean = nullptr;
        MonitorElement *h_stddev = nullptr;

        Profile(DQMStore::IBooker &iBooker, MonitorElement *_h);

        void Fill(double x, double y);
        void Write() const;
    };

    struct SectorData
    {

    };

    // ------------ member data ------------
    edm::EDGetTokenT<CTPPSLocalTrackLiteCollection> tracksToken_;
};

// -------------------------------- Stat methods --------------------------------
PPSAlignmentWorker::Stat::Stat() {}
PPSAlignmentWorker::Stat::Stat(unsigned int dim_)
{
    init(dim_);
}

void PPSAlignmentWorker::Stat::init(unsigned int dim_)
{
    dim = dim_;
    S1 = 0.;

    Sv.resize(dim);
    Svv.resize(dim);
    Svvv.resize(dim);
    Svvvv.resize(dim);

    std::vector<double> tmp(dim);
    Sxy.resize(dim, tmp);
    Sxxy.resize(dim, tmp);
    Sxyy.resize(dim, tmp);
    Sxxyy.resize(dim, tmp);
}

void PPSAlignmentWorker::Stat::setLabels(const std::vector<std::string> &labels_)
{
    labels.resize(dim);
    for (unsigned int i = 0; i < dim; i++)
        labels[i] = labels_[i];
}

template<class T>
void PPSAlignmentWorker::Stat::fill(const T &v)
{
    S1 += 1.;
    for (unsigned int i = 0; i < dim; i++)
    {
        Sv[i] += v[i];
        Svv[i] += v[i] * v[i];
        Svvv[i] += v[i] * v[i] * v[i];
        Svvvv[i] += v[i] * v[i] * v[i] * v[i];

        for (unsigned int j = 0; j < dim; j++)
        {
            Sxy[i][j] += v[i] * v[j];
            Sxxy[i][j] += v[i] * v[i] * v[j];
            Sxyy[i][j] += v[i] * v[j] * v[j];
            Sxxyy[i][j] += v[i] * v[i] * v[j] * v[j];
        }
    }
}

void PPSAlignmentWorker::Stat::fill(double v1, double v2, double v3, double v4, double v5)
{
    std::vector<double> v = {v1, v2, v3, v4, v5};
    fill(v);
}

std::string PPSAlignmentWorker::Stat::qLabel(unsigned int i) const
{
    if (labels.empty())
    {
        char buf[10];
        sprintf(buf, "qu.%3i", i);
        return buf;
    }
    else
    {
        return labels[i];
    }
}

// --- 1D getters ---
double PPSAlignmentWorker::Stat::entries() const
{
    return S1;
}

double PPSAlignmentWorker::Stat::mean(unsigned int i) const
{
    return (S1 > 0.) ? Sv[i] / S1 : 0.;
}

double PPSAlignmentWorker::Stat::stdDev(unsigned int i) const
{
    double v = (Svv[i] - Sv[i] * Sv[i] / S1) / (S1- 1.);
    return (v > 0.) ? std::sqrt(v) : 0.;
}


double PPSAlignmentWorker::Stat::meanUnc(unsigned int i) const
{
    return (S1 > 0.) ? stdDev(i) / sqrt(S1) : 0.;
}

double PPSAlignmentWorker::Stat::stdDevUnc(unsigned int i) const
{
    double mu = mean(i);
    double s = stdDev(i);
    double v = s*s;

    double sum = Svvvv[i] - 4. * mu * Svvv[i] + 6. * mu*mu * Svv[i] - 4. * mu*mu*mu * Sv[i] + mu*mu*mu*mu * S1;
    double E4 = (S1 > 1.) ? sum / (S1 - 1.) : 0.;

    double v_var  = (S1 > 3.) ? (E4 - (S1 - 3.) / (S1 - 1.) * v*v) / S1 : 0.;
    double s_var = v_var / 4. / v;
    return (s_var > 0.) ? std::sqrt(s_var) : 0.;
}

// approximation of GetStdDevUnc valid for Gaussian distributions
double PPSAlignmentWorker::Stat::stdDevUncGauss(unsigned int i) const
{
    double s = stdDev(i);
    return (S1 > 0.) ? s / std::sqrt(2. * S1) : 0.;
}

// --- 2D getters ---
double PPSAlignmentWorker::Stat::covariance(unsigned int i, unsigned int j) const
{
    return (S1 > 1.) ? (Sxy[i][j] - Sv[i] * Sv[j] / S1) / (S1 - 1.) : 0.;
}

double PPSAlignmentWorker::Stat::correlation(unsigned int i, unsigned int j) const
{
    double C = covariance(i, j);
    double den = stdDev(i) * stdDev(j);
    return (den > 0.) ? C / den : 0.;
}

double PPSAlignmentWorker::Stat::covarianceUnc(unsigned int i, unsigned int j) const
{
    double mx = mean(i);
    double my = mean(j);
    double sx = stdDev(i);
    double sy = stdDev(j);
    double C = covariance(i, j);

    double sum =
        Sxxyy[i][j] 
        -2.*Sxyy[i][j]*mx - 2.*Sxxy[i][j]*my
        + 4.*Sxy[i][j]*mx*my
        + Svv[i]*my*my + Svv[j]*mx*mx
        - 2.*Sv[i]*mx*my*my - 2.*Sv[j]*mx*mx*my
        + mx*mx*my*my;
    double D = (S1 > 1.) ? sum / (S1 - 1.) : 0.;

    double C_var = (S1 > 2.) ? (D + sx*sx*sy*sy/(S1 - 1.) - (S1-2.)/(S1-1.)*C*C) / S1 : 0.;
    return (C_var > 0.) ? std::sqrt(C_var) : 0.;
}

double PPSAlignmentWorker::Stat::correlationUnc(unsigned int i, unsigned int j) const
{
    // WARNING: the calculation below assumes no correlation between C, si_i and si_j, which
    // might not be correct - in that case it gives an upper bound for the uncertainty
    double C = covariance(i, j), C_unc = covarianceUnc(i, j);
    double si_i = stdDev(i), si_i_unc = stdDevUnc(i);
    double si_j = stdDev(j), si_j_unc = stdDevUnc(j);
    double rho = C / (si_i * si_j);
    double sum =
        (C != 0. && si_i != 0. && si_j != 0.) ? std::pow(C_unc / C, 2.) + std::pow(si_i_unc / si_i, 2.) + std::pow(si_j_unc / si_j, 2.) : 0.;
    double rho_unc = std::fabs(rho) * std::sqrt(sum);
    return rho_unc;
}

TMatrixDSym PPSAlignmentWorker::Stat::covarianceMatrix() const
{
    TMatrixDSym m(dim);
    for (unsigned int i = 0; i < dim; i++)
    {
        for (unsigned int j = 0; j < dim; j++)
            m(i, j) = covariance(i, j);
    }
    return m;
}

// --- print methods ---
void PPSAlignmentWorker::Stat::printStat() const
{
    printf("entries: %.3E\n", S1);
}

void PPSAlignmentWorker::Stat::printMeanAndStdDev() const
{
    for (unsigned int i = 0; i < dim; i++)
    {
        double mu = mean(i);
        double mu_unc = meanUnc(i);
        double s = stdDev(i);
        double s_unc = stdDevUnc(i);
        printf("%s: mean %+.3E +- %.3E, std. dev. = %.3E +- %.3E\n", qLabel(i).c_str(), mu, mu_unc, s, s_unc);
    }
}

void PPSAlignmentWorker::Stat::printCovariance() const
{
    printf("      ");
    for (unsigned int i = 0; i < dim; i++)
        printf("   %6s", qLabel(i).c_str());

    printf("\n");

    for (unsigned int i = 0; i < dim; i++)
    {
        printf("%6s", qLabel(i).c_str());
        for (unsigned int j = 0; j < dim; j++)
            printf("   %+.3f", covariance(i, j));
        printf("\n");
    }
}

void PPSAlignmentWorker::Stat::printCorrelation() const
{
    printf("      ");
    for (unsigned int i = 0; i < dim; i++)
        printf("   %6s", qLabel(i).c_str());

    printf("\n");

    for (unsigned int i = 0; i < dim; i++)
    {
        printf("%6s", qLabel(i).c_str());
        for (unsigned int j = 0; j < dim; j++)
            printf("   %+.3f", correlation(i, j));

        printf("\n");
    }
}

// -------------------------------- PPSAlignmentWorker methods --------------------------------

PPSAlignmentWorker::PPSAlignmentWorker(const edm::ParameterSet &iConfig) 
    : tracksToken_(consumes<CTPPSLocalTrackLiteCollection>(iConfig.getParameter<edm::InputTag>("tagTracks")))
{}

// ------------ method called for each event  ------------
void PPSAlignmentWorker::analyze(const edm::Event &iEvent, const edm::EventSetup &iSetup)
{
    edm::ESHandle<PPSAlignmentConfig> config;
    iSetup.get<PPSAlignmentConfigRcd>().get(config);

    edm::Handle<CTPPSLocalTrackLiteCollection> tracks;
    iEvent.getByToken(tracksToken_, tracks);
}


void PPSAlignmentWorker::bookHistograms(DQMStore::IBooker &ibook, edm::Run const &run, edm::EventSetup const &iSetup)
{

}

DEFINE_FWK_MODULE(PPSAlignmentWorker);