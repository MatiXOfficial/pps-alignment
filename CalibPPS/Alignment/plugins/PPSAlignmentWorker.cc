/****************************************************************************
 *
 *  CalibPPS/Alignment/plugins/PPSAlignmentWorker.cc
 *
 *  Description : PPS Alignment DQM worker
 *
 * Authors:
 *  - Jan Kašpar
 *  - Mateusz Kocot
 *
 ****************************************************************************/

#include "DQMServices/Core/interface/DQMOneEDAnalyzer.h"
#include "DQMServices/Core/interface/DQMStore.h"
#include "DQMServices/Core/interface/MonitorElement.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/CTPPSDetId/interface/TotemRPDetId.h"
#include "DataFormats/CTPPSReco/interface/CTPPSLocalTrackLite.h"
#include "DataFormats/CTPPSReco/interface/CTPPSLocalTrackLiteFwd.h"

#include "CondFormats/PPSObjects/interface/PPSAlignmentConfig.h"
#include "CondFormats/DataRecord/interface/PPSAlignmentConfigRcd.h"

#include <vector>
#include <string>
#include <cmath>

#include "TVectorD.h"
#include "TH2D.h"
#include "TGraph.h"

//----------------------------------------------------------------------------------------------------

class PPSAlignmentWorker : public DQMOneEDAnalyzer<>
{
public:
    PPSAlignmentWorker(const edm::ParameterSet &);
    ~PPSAlignmentWorker() override {};

private:
    void bookHistograms(DQMStore::IBooker &, edm::Run const &, edm::EventSetup const &) override;
    void analyze(const edm::Event &, const edm::EventSetup &) override;
    void dqmEndRun(edm::Run const &, edm::EventSetup const &) override;

    // ------------ structures ------------
    struct Stat
    {
        unsigned int dim;
        double S1;
        std::vector<double> Sv, Svv, Svvv, Svvvv;
        std::vector< std::vector<double> > Sxy, Sxxy, Sxyy, Sxxyy;
        std::vector<std::string> labels;

        Stat();
        Stat(unsigned int _dim);

        void init(unsigned int _dim = 1);
        void setLabels(const std::vector<std::string> &_labels);

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
        TH2D *h;
        std::vector<Stat> st;

        MonitorElement *h_entries;
        MonitorElement *h_mean;
        MonitorElement *h_stddev;

        Profile();
        Profile(DQMStore::IBooker &iBooker, TH2D *_h);
        void fill(double x, double y);
        void fillEndRun() const;
    };

    struct SectorData
    {
        std::string name;

        unsigned int rpIdUp, rpIdDw;

        SectorConfig scfg;

        // hit distributions
        std::map<unsigned int, MonitorElement*> m_h1_x_bef_sel;

        std::map<unsigned int, MonitorElement*> m_h2_y_vs_x_bef_sel;

        std::map<unsigned int, MonitorElement*> m_h2_y_vs_x_mlt_sel;
        
        std::map<unsigned int, MonitorElement*> m_h2_y_vs_x_aft_sel;
        std::map<unsigned int, TGraph*> m_g_y_vs_x_aft_sel;

        // cut plots
        MonitorElement *h_q_cut_h_bef, *h_q_cut_h_aft;
        MonitorElement *h2_cut_h_bef, *h2_cut_h_aft;
        MonitorElement *p_cut_h_aft;

        MonitorElement *h_q_cut_v_bef, *h_q_cut_v_aft;
        MonitorElement *h2_cut_v_bef, *h2_cut_v_aft;
        MonitorElement *p_cut_v_aft;

        // profiles
        std::map<unsigned int, Profile> m_p_y_vs_x_aft_sel;

        // near-far plots
        MonitorElement *p_x_diffFN_vs_x_N;
        MonitorElement *p_y_diffFN_vs_y_N;
        MonitorElement *p_y_diffFN_vs_y_F;

        struct SlicePlots
        {
            MonitorElement *h_y;
            MonitorElement *h2_y_diffFN_vs_y;
            MonitorElement *p_y_diffFN_vs_y;

            SlicePlots();
            SlicePlots(DQMStore::IBooker &iBooker);
        };

        std::map<unsigned int, SlicePlots> x_slice_plots_N, x_slice_plots_F;

        void init(DQMStore::IBooker &iBooker, edm::EventSetup const &iSetup, std::string _name, 
                  unsigned int _rpIdUp, unsigned int _rpIdDw, const SectorConfig &_scfg, const std::string &folder);

        unsigned int process(const edm::EventSetup &iSetup, const CTPPSLocalTrackLiteCollection &tracks);
        void fillEndRun() const;
    };

    // ------------ member data ------------
    edm::EDGetTokenT<CTPPSLocalTrackLiteCollection> tracksToken_;

    SectorData sectorData45;
    SectorData sectorData56;

    unsigned long eventCount = 0;
    unsigned long eventSelCount45 = 0;
    unsigned long eventSelCount56 = 0;

    std::string folder_;
};

// -------------------------------- Stat methods --------------------------------

PPSAlignmentWorker::Stat::Stat() {}
PPSAlignmentWorker::Stat::Stat(unsigned int _dim)
{
    init(_dim);
}

void PPSAlignmentWorker::Stat::init(unsigned int _dim)
{
    dim = _dim;
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

void PPSAlignmentWorker::Stat::setLabels(const std::vector<std::string> &_labels)
{
    labels.resize(dim);
    for (unsigned int i = 0; i < dim; i++)
        labels[i] = _labels[i];
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

// -------------------------------- Profile methods --------------------------------
PPSAlignmentWorker::Profile::Profile()  {}

PPSAlignmentWorker::Profile::Profile(DQMStore::IBooker &iBooker, TH2D *_h)
    : h(_h)
{
    st.resize(h->GetNbinsX(), Stat(1));

    const int bins = h->GetXaxis()->GetNbins();
    const double xMin = h->GetXaxis()->GetXmin();
    const double xMax = h->GetXaxis()->GetXmax();

    h_entries = iBooker.book1D("h_entries", ";x", bins, xMin, xMax);
    h_mean = iBooker.book1D("h_mean", ";x", bins, xMin, xMax);
    h_stddev = iBooker.book1D("h_sttdev", ";x", bins, xMin, xMax);
}

void PPSAlignmentWorker::Profile::fill(double x, double y)
{
    int bi = h->GetXaxis()->FindBin(x);

    if (bi < 1 || bi > h->GetNbinsX())
        return;

    int vi = bi - 1;
    st[vi].fill(y);
}

void PPSAlignmentWorker::Profile::fillEndRun() const
{
    const int bins = h->GetXaxis()->GetNbins();

    for (int bi = 1; bi <= bins; bi++)
    {
        int vi = bi - 1;

        h_entries->setBinContent(bi, st[vi].entries());

        h_mean->setBinContent(bi, st[vi].mean(0));
        h_mean->setBinError(bi, st[vi].meanUnc(0));

        h_stddev->setBinContent(bi, st[vi].stdDev(0));
        h_stddev->setBinError(bi, st[vi].stdDevUncGauss(0));
    }
}

// -------------------------------- SectorData and SlicePlots methods --------------------------------

PPSAlignmentWorker::SectorData::SlicePlots::SlicePlots() {}

PPSAlignmentWorker::SectorData::SlicePlots::SlicePlots(DQMStore::IBooker &iBooker)
{
    h_y = iBooker.book1D("h_y", ";y", 100, -10., 10.);
    h2_y_diffFN_vs_y = iBooker.book2DD("h2_y_diffFN_vs_y", ";y;x_{F} - y_{N}", 100, -10., 10., 100, -2., 2.);
    auto tmp_p_y_diffFN_vs_y = new TProfile("", ";y;x_{F} - y_{N}", 100, -10., 10.);
    p_y_diffFN_vs_y = iBooker.bookProfile("p_y_diffFN_vs_y", tmp_p_y_diffFN_vs_y);
}

void PPSAlignmentWorker::SectorData::init(DQMStore::IBooker &iBooker, edm::EventSetup const &iSetup, std::string _name, 
                                          unsigned int _rpIdUp, unsigned int _rpIdDw, const SectorConfig &_scfg, const std::string &folder)
{
    edm::ESHandle<PPSAlignmentConfig> cfg;
    iSetup.get<PPSAlignmentConfigRcd>().get(cfg);

    name = _name;
    rpIdUp = _rpIdUp;
    rpIdDw = _rpIdDw;
    scfg = _scfg;

    // binning
    const double bin_size_x = 142.3314E-3; // mm
	const unsigned int n_bins_x = 210;

	const double pixel_x_offset = (cfg->aligned()) ? 0. : 40.;

	const double x_min_pix = pixel_x_offset, x_max_pix = pixel_x_offset + n_bins_x * bin_size_x;
	const double x_min_str = 0., x_max_str = n_bins_x * bin_size_x;

	const unsigned int n_bins_y = 400;
	const double y_min = -20., y_max = +20.;

    // hit distributions
    iBooker.setCurrentFolder(folder + "/worker/" + _name + "/before selection/" + cfg->rpTags()[rpIdUp]);
    m_h1_x_bef_sel[rpIdUp] = iBooker.book1DD("h_x", ";x", 10 * n_bins_x, x_min_str, x_max_str);
	m_h2_y_vs_x_bef_sel[rpIdUp] = iBooker.book2DD("h2_y_vs_x", ";x;y", n_bins_x, x_min_str, x_max_str, n_bins_y, y_min, y_max);
    iBooker.setCurrentFolder(folder + "/worker/" + _name + "/before selection/" + cfg->rpTags()[rpIdDw]);
	m_h2_y_vs_x_bef_sel[rpIdDw] = iBooker.book2DD("h2_y_vs_x", ";x;y", n_bins_x, x_min_pix, x_max_pix, n_bins_y, y_min, y_max);
    m_h1_x_bef_sel[rpIdDw] = iBooker.book1DD("h_x", ";x", 10 * n_bins_x, x_min_pix, x_max_pix);

    iBooker.setCurrentFolder(folder + "/worker/" + _name + "/multiplicity selection/" + cfg->rpTags()[rpIdUp]);
	m_h2_y_vs_x_mlt_sel[rpIdUp] = iBooker.book2DD("h2_y_vs_x", ";x;y", n_bins_x, x_min_str, x_max_str, n_bins_y, y_min, y_max);
    iBooker.setCurrentFolder(folder + "/worker/" + _name + "/multiplicity selection/" + cfg->rpTags()[rpIdDw]);
	m_h2_y_vs_x_mlt_sel[rpIdDw] = iBooker.book2DD("h2_y_vs_x", ";x;y", n_bins_x, x_min_pix, x_max_pix, n_bins_y, y_min, y_max);

    iBooker.setCurrentFolder(folder + "/worker/" + _name + "/after selection/" + cfg->rpTags()[rpIdUp]);
	m_h2_y_vs_x_aft_sel[rpIdUp] = iBooker.book2DD("h2_y_vs_x", ";x;y", n_bins_x, x_min_str, x_max_str, n_bins_y, y_min, y_max);
    iBooker.setCurrentFolder(folder + "/worker/" + _name + "/after selection/" + cfg->rpTags()[rpIdDw]);
	m_h2_y_vs_x_aft_sel[rpIdDw] = iBooker.book2DD("h2_y_vs_x", ";x;y", n_bins_x, x_min_pix, x_max_pix, n_bins_y, y_min, y_max);

	m_g_y_vs_x_aft_sel[rpIdUp] = new TGraph();
	m_g_y_vs_x_aft_sel[rpIdDw] = new TGraph();

    // cut plots
    iBooker.setCurrentFolder(folder + "/worker/" + _name + "/cuts/cut_h");
    h_q_cut_h_bef = iBooker.book1D("h_q_cut_h_bef", ";cq_h", 400, -2., 2.);
	h_q_cut_h_aft = iBooker.book1D("h_q_cut_h_aft", ";cq_h", 400, -2., 2.);
	h2_cut_h_bef = iBooker.book2DD("h2_cut_h_bef", ";x_up;x_dw", n_bins_x, x_min_str, x_max_str, n_bins_x, x_min_pix, x_max_pix);
	h2_cut_h_aft = iBooker.book2DD("h2_cut_h_aft", ";x_up;x_dw", n_bins_x, x_min_str, x_max_str, n_bins_x, x_min_pix, x_max_pix);
	auto tmp_p_cut_h_aft = new TProfile("", ";x_up;mean of x_dw", n_bins_x, x_min_str, x_max_str);
    p_cut_h_aft = iBooker.bookProfile("p_cut_h_aft", tmp_p_cut_h_aft);

    iBooker.setCurrentFolder(folder + "/worker/" + _name + "/cuts/cut_v");
	h_q_cut_v_bef = iBooker.book1D("h_q_cut_v_bef", ";cq_v", 400, -2., 2.);
	h_q_cut_v_aft = iBooker.book1D("h_q_cut_v_aft", ";cq_v", 400, -2., 2.);
	h2_cut_v_bef = iBooker.book2DD("h2_cut_v_bef", ";y_up;y_dw", n_bins_y, y_min, y_max, n_bins_y, y_min, y_max);
	h2_cut_v_aft = iBooker.book2DD("h2_cut_v_aft", ";y_up;y_dw", n_bins_y, y_min, y_max, n_bins_y, y_min, y_max);
	auto tmp_p_cut_v_aft = new TProfile("", ";y_up;mean of y_dw", n_bins_y, y_min, y_max);
    p_cut_v_aft = iBooker.bookProfile("p_cut_v_aft", tmp_p_cut_v_aft);

    // profiles
    iBooker.setCurrentFolder(folder + "/worker/" + _name + "/profiles/" + cfg->rpTags()[rpIdUp]);
    m_p_y_vs_x_aft_sel.insert({rpIdUp, Profile(iBooker, m_h2_y_vs_x_aft_sel[rpIdUp]->getTH2D())});
    iBooker.setCurrentFolder(folder + "/worker/" + _name + "/profiles/" + cfg->rpTags()[rpIdDw]);
    m_p_y_vs_x_aft_sel.insert({rpIdDw, Profile(iBooker, m_h2_y_vs_x_aft_sel[rpIdDw]->getTH2D())});

    // near-far plots
    iBooker.setCurrentFolder(folder + "/worker/" + _name + "/near_far");
	auto tmp_p_x_diffFN_vs_x_N = new TProfile("", ";x_{N};x_{F} - x_{N}", 100, 0., 20.);
    p_x_diffFN_vs_x_N = iBooker.bookProfile("p_x_diffFN_vs_x_N", tmp_p_x_diffFN_vs_x_N);
	auto tmp_p_y_diffFN_vs_y_N = new TProfile("", ";y_{N};y_{F} - y_{N}", 200, -10., 10.);
    p_y_diffFN_vs_y_N = iBooker.bookProfile("p_y_diffFN_vs_y_N", tmp_p_y_diffFN_vs_y_N);
	auto tmp_p_y_diffFN_vs_y_F = new TProfile("", ";y_{F};y_{F} - y_{N}", 200, -10., 10.);
    p_y_diffFN_vs_y_F = iBooker.bookProfile("p_y_diffFN_vs_y_F", tmp_p_y_diffFN_vs_y_F);

    for (int i = 0; i < scfg.nr_x_slice_n; i++)
    {
        const double xMin = scfg.nr_x_slice_min + i * scfg.nr_x_slice_w;
        const double xMax = scfg.nr_x_slice_min + (i + 1) * scfg.nr_x_slice_w;

        char buf[100];
        sprintf(buf, "%.1f-%.1f", xMin, xMax);

        iBooker.setCurrentFolder(folder + "/worker/" + _name + "/near_far/x slices, N/" + buf);
        x_slice_plots_N.insert({i, SlicePlots(iBooker)});
    }

    for (int i = 0; i < scfg.fr_x_slice_n; i++)
    {
        const double xMin = scfg.fr_x_slice_min + i * scfg.fr_x_slice_w;
        const double xMax = scfg.fr_x_slice_min + (i + 1) * scfg.fr_x_slice_w;

        char buf[100];
        sprintf(buf, "%.1f-%.1f", xMin, xMax);

        iBooker.setCurrentFolder(folder + "/worker/" + _name + "/near_far/x slices, F/" + buf);
        x_slice_plots_F.insert({i, SlicePlots(iBooker)});
    }
}

unsigned int PPSAlignmentWorker::SectorData::process(const edm::EventSetup &iSetup, const CTPPSLocalTrackLiteCollection &tracks)
{
    edm::ESHandle<PPSAlignmentConfig> cfg;
    iSetup.get<PPSAlignmentConfigRcd>().get(cfg);

    CTPPSLocalTrackLiteCollection tracksUp, tracksDw;

    for (const auto &tr : tracks)
	{
		CTPPSDetId rpId(tr.rpId());
		unsigned int rpDecId = rpId.arm()*100 + rpId.station()*10 + rpId.rp();

        if (rpDecId != rpIdUp && rpDecId != rpIdDw)
			continue;

		double x = tr.x();
		double y = tr.y();

        // apply alignment corrections
		x += cfg->alignmentCorrectionsX()[rpDecId];
		y += cfg->alignmentCorrectionsY()[rpDecId];

        // re-build track object
        CTPPSLocalTrackLite trCorr(tr.rpId(), x, 0., y, 0.,
        tr.tx(), tr.txUnc(), tr.ty(), tr.tyUnc(),
        tr.chiSquaredOverNDF(), tr.pixelTrackRecoInfo(), tr.numberOfPointsUsedForFit(),
        tr.time(), tr.timeUnc());

        // store corrected track into the right collection
        if (rpDecId == rpIdUp)
			tracksUp.push_back(std::move(trCorr));
		if (rpDecId == rpIdDw)
			tracksDw.push_back(std::move(trCorr));
    }

    // update plots before selection
    for (const auto &tr : tracksUp)
	{
		m_h1_x_bef_sel[rpIdUp]->Fill(tr.x());
		m_h2_y_vs_x_bef_sel[rpIdUp]->Fill(tr.x(), tr.y());
	}

    for (const auto &tr : tracksDw)
	{
		m_h1_x_bef_sel[rpIdDw]->Fill(tr.x());
		m_h2_y_vs_x_bef_sel[rpIdDw]->Fill(tr.x(), tr.y());
	}

    // skip crowded events
    if (tracksUp.size() > 2)
		return 0;

	if (tracksDw.size() > 2)
		return 0;

    // update plots with multiplicity selection
	for (const auto &tr : tracksUp)
		m_h2_y_vs_x_mlt_sel[rpIdUp]->Fill(tr.x(), tr.y());

	for (const auto &tr : tracksDw)
		m_h2_y_vs_x_mlt_sel[rpIdDw]->Fill(tr.x(), tr.y());

    // do the selection
    unsigned int pairsSelected = 0;

    for (const auto &trUp : tracksUp)
    {
        for (const auto &trDw : tracksDw)
        {
            h2_cut_h_bef->Fill(trUp.x(), trDw.x());
			h2_cut_v_bef->Fill(trUp.y(), trDw.y());

            const double cq_h = trDw.x() + scfg.cut_h_a * trUp.x() + scfg.cut_h_c;
			h_q_cut_h_bef->Fill(cq_h);
			const bool cv_h = (std::fabs(cq_h) < cfg->n_si() * scfg.cut_h_si);

			const double cq_v = trDw.y() + scfg.cut_v_a * trUp.y() + scfg.cut_v_c;
			h_q_cut_v_bef->Fill(cq_v);
			const bool cv_v = (std::fabs(cq_v) < cfg->n_si() * scfg.cut_v_si);

			bool cutsPassed = true;
			if (scfg.cut_h_apply)
				cutsPassed &= cv_h;
			if (scfg.cut_v_apply)
				cutsPassed &= cv_v;

            if (cutsPassed)
            {
                pairsSelected++;

				h_q_cut_h_aft->Fill(cq_h);
				h_q_cut_v_aft->Fill(cq_v);

				h2_cut_h_aft->Fill(trUp.x(), trDw.x());
				h2_cut_v_aft->Fill(trUp.y(), trDw.y());

				p_cut_h_aft->Fill(trUp.x(), trDw.x());
				p_cut_v_aft->Fill(trUp.y(), trDw.y());

				m_h2_y_vs_x_aft_sel[rpIdUp]->Fill(trUp.x(), trUp.y());
				m_h2_y_vs_x_aft_sel[rpIdDw]->Fill(trDw.x(), trDw.y());

				int idx = m_g_y_vs_x_aft_sel[rpIdUp]->GetN();
				m_g_y_vs_x_aft_sel[rpIdUp]->SetPoint(idx, trUp.x(), trUp.y());
				m_g_y_vs_x_aft_sel[rpIdDw]->SetPoint(idx, trDw.x(), trDw.y());

				m_p_y_vs_x_aft_sel[rpIdUp].fill(trUp.x(), trUp.y());
				m_p_y_vs_x_aft_sel[rpIdDw].fill(trDw.x(), trDw.y());

				p_x_diffFN_vs_x_N->Fill(trUp.x(), trDw.x() - trUp.x());

                const auto &range = cfg->alignment_y_ranges()[rpIdUp];          // probably redundant    
				if (trUp.x() > range.x_min && trUp.x() < range.x_max)           //
				{                                                               //
					p_y_diffFN_vs_y_N->Fill(trUp.y(), trDw.y() - trUp.y());     //
					p_y_diffFN_vs_y_F->Fill(trDw.y(), trDw.y() - trUp.y());     //
				}                                                               //

				idx = (trUp.x() - scfg.nr_x_slice_min) / scfg.nr_x_slice_w;
				if (idx >= 0 && idx < scfg.nr_x_slice_n)
				{
					x_slice_plots_N[idx].h_y->Fill(trUp.y());
					x_slice_plots_N[idx].h2_y_diffFN_vs_y->Fill(trUp.y(), trDw.y() - trUp.y());
					x_slice_plots_N[idx].p_y_diffFN_vs_y->Fill(trUp.y(), trDw.y() - trUp.y());
				}

				idx = (trDw.x() - scfg.fr_x_slice_min) / scfg.fr_x_slice_w;
				if (idx >= 0 && idx < scfg.fr_x_slice_n)
				{
					x_slice_plots_F[idx].h_y->Fill(trDw.y());
					x_slice_plots_F[idx].h2_y_diffFN_vs_y->Fill(trDw.y(), trDw.y() - trUp.y());
					x_slice_plots_F[idx].p_y_diffFN_vs_y->Fill(trDw.y(), trDw.y() - trUp.y());
				}
            }
        }
    }

    return pairsSelected;
}

void PPSAlignmentWorker::SectorData::fillEndRun() const
{
    for (const auto &p : m_p_y_vs_x_aft_sel)
    {
        p.second.fillEndRun();
    }
}

// -------------------------------- PPSAlignmentWorker methods --------------------------------

PPSAlignmentWorker::PPSAlignmentWorker(const edm::ParameterSet &iConfig) 
    : tracksToken_(consumes<CTPPSLocalTrackLiteCollection>(iConfig.getParameter<edm::InputTag>("tagTracks"))),
      folder_(iConfig.getParameter<std::string>("folder"))
{}

// ------------ method called for each event  ------------
void PPSAlignmentWorker::analyze(const edm::Event &iEvent, const edm::EventSetup &iSetup)
{
    const auto &tracks = iEvent.get(tracksToken_);

    sectorData45.process(iSetup, tracks);
    sectorData56.process(iSetup, tracks);
}


void PPSAlignmentWorker::bookHistograms(DQMStore::IBooker &iBooker, edm::Run const &, edm::EventSetup const &iSetup)
{
    edm::ESHandle<PPSAlignmentConfig> cfg;
    iSetup.get<PPSAlignmentConfigRcd>().get(cfg);

    sectorData45.init(iBooker, iSetup, "sector 45", 3, 23, cfg->sectorConfig45(), folder_);
    sectorData56.init(iBooker, iSetup, "sector 56", 103, 123, cfg->sectorConfig56(), folder_);
}

void PPSAlignmentWorker::dqmEndRun(edm::Run const &, edm::EventSetup const &)
{
    sectorData45.fillEndRun();
    sectorData56.fillEndRun();
}

DEFINE_FWK_MODULE(PPSAlignmentWorker);