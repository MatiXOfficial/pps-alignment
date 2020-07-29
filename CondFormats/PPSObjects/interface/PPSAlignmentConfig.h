/****************************************************************************
 *
 *  CondFormats/PPSObjects/interface/PPSAlignmentConfig.h
 *
 *  Description : Alignment parameters
 *
 * Authors:
 *  - Jan Kašpar
 *  - Mateusz Kocot
 *
 ****************************************************************************/

#ifndef CondFormats_PPSObjects_PPSAlignmentConfig_h
#define CondFormats_PPSObjects_PPSAlignmentConfig_h

#include "CondFormats/Serialization/interface/Serializable.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include <vector>
#include <string>
#include <map>
#include <cmath>
#include <iomanip>
#include <iostream>

//---------------------------------------------------------------------------------------------

struct SelectionRange
{
    double x_min;
    double x_max;

    SelectionRange(double x_min = 0., double x_max = 0.) : x_min(x_min), x_max(x_max) {}

    COND_SERIALIZABLE;
};

//---------------------------------------------------------------------------------------------

struct SectorConfig
{
    bool cut_h_apply;
    double cut_h_a, cut_h_c, cut_h_si;

    bool cut_v_apply;
    double cut_v_a, cut_v_c, cut_v_si;

    int nr_x_slice_n;
	double nr_x_slice_min, nr_x_slice_w;
	int fr_x_slice_n;
	double fr_x_slice_min, fr_x_slice_w;

    COND_SERIALIZABLE;
};

//---------------------------------------------------------------------------------------------

class PPSAlignmentConfig
{
public:
    // Getters
    unsigned int fill() const { return fill_; }
    unsigned int xangle() const { return xangle_; }
    double beta() const { return beta_; }
    std::string dataset() const { return dataset_; }

	std::map<unsigned int, std::string> rpTags() const { return rpTags_; }

	std::vector<std::string> inputFiles() const { return inputFiles_; }

	std::map<unsigned int, double> alignmentCorrectionsX() const { return alignmentCorrectionsX_; }
    std::map<unsigned int, double> alignmentCorrectionsY() const { return alignmentCorrectionsY_; }

	bool aligned() const { return aligned_; }

	double n_si() const { return n_si_; }

	SectorConfig sectorConfig45() const { return sectorConfig45_; }
    SectorConfig sectorConfig56() const { return sectorConfig56_; }

	std::vector<std::string> matchingReferenceDatasets() const { return matchingReferenceDatasets_; }
	std::map<unsigned int, SelectionRange> matchingShiftRanges() const { return matchingShiftRanges_; }

	std::map<unsigned int, SelectionRange> alignment_x_meth_o_ranges() const { return alignment_x_meth_o_ranges_; }
	std::map<unsigned int, SelectionRange> alignment_x_relative_ranges() const { return alignment_x_relative_ranges_; }

	std::map<unsigned int, SelectionRange> alignment_y_ranges() const { return alignment_y_ranges_; }

    // Setters 
    void setFill(unsigned int fill) { fill_ = fill; }
    void setXangle(unsigned int xangle) { xangle_ = xangle; }
    void setBeta(double beta) { beta_ = beta; }
    void setDataset(std::string &dataset) { dataset_ = dataset; }

    void setRpTags(std::map<unsigned int, std::string> &rpTags) { rpTags_ = rpTags; }

    void setInputFiles(std::vector<std::string> &inputFiles) { inputFiles_ = inputFiles; }

    void setAlignmentCorrectionsX(std::map<unsigned int, double> &alignmentCorrectionsX) 
    { 
        alignmentCorrectionsX_ = alignmentCorrectionsX; 
    }
    void setAlignmentCorrectionsY(std::map<unsigned int, double> &alignmentCorrectionsY) 
    { 
        alignmentCorrectionsY_ = alignmentCorrectionsY; 
    }

    void setAligned(bool aligned) { aligned_ = aligned; }

    void setN_si(double n_si) { n_si_ = n_si; }

    void setSectorConfig45(SectorConfig &sectorConfig45) { sectorConfig45_ = sectorConfig45; }
    void setSectorConfig56(SectorConfig &sectorConfig56) { sectorConfig56_ = sectorConfig56; }

    void setMatchingReferenceDatasets(std::vector<std::string> &matchingReferenceDatasets)
    {
        matchingReferenceDatasets_ = matchingReferenceDatasets;
    }
    void setMatchingShiftRanges(std::map<unsigned int, SelectionRange> &matchingShiftRanges)
    {
        matchingShiftRanges_ = matchingShiftRanges;
    }

    void setAlignment_x_meth_o_ranges(std::map<unsigned int, SelectionRange> &alignment_x_meth_o_ranges)
    {
        alignment_x_meth_o_ranges_ = alignment_x_meth_o_ranges;
    }
    void setAlignment_x_relative_ranges(std::map<unsigned int, SelectionRange> &alignment_x_relative_ranges)
    {
        alignment_x_relative_ranges_ = alignment_x_relative_ranges;   
    }

    void setAlignment_y_ranges(std::map<unsigned int, SelectionRange> &alignment_y_ranges)
    {
        alignment_y_ranges_ = alignment_y_ranges;
    }

    // << operator
    friend std::ostream &operator<<(std::ostream &os, PPSAlignmentConfig c);

private:
    unsigned int fill_;
	unsigned int xangle_;
	double beta_;
	std::string dataset_;

	std::map<unsigned int, std::string> rpTags_;

	std::vector<std::string> inputFiles_;

	std::map<unsigned int, double> alignmentCorrectionsX_, alignmentCorrectionsY_;

	bool aligned_;

	double n_si_;

	SectorConfig sectorConfig45_, sectorConfig56_;

	std::vector<std::string> matchingReferenceDatasets_;
	std::map<unsigned int, SelectionRange> matchingShiftRanges_;

	std::map<unsigned int, SelectionRange> alignment_x_meth_o_ranges_;
	std::map<unsigned int, SelectionRange> alignment_x_relative_ranges_;

	std::map<unsigned int, SelectionRange> alignment_y_ranges_;

    COND_SERIALIZABLE;
};

//---------------------------------------------------------------------------------------------

std::ostream &operator<<(std::ostream &os, SectorConfig &sc)
{
    os << std::fixed << std::setprecision(3);
    os << "    cut_h: apply = " << sc.cut_h_apply << ", a = " << sc.cut_h_a << ", c = " 
        << sc.cut_h_c << ", si = " << sc.cut_h_si << "\n";
    os << "    cut_v: apply = " << sc.cut_v_apply << ", a = " << sc.cut_v_a << ", c = " 
        << sc.cut_v_c << ", si = " << sc.cut_v_si << "\n";
    os << std::setprecision(2);
    os << "    x slices, nr: min = " << sc.nr_x_slice_min << ", w = " << sc.nr_x_slice_w 
        << ", n = " << sc.nr_x_slice_n << "\n";
    os << "    x slices, fr: min = " << sc.fr_x_slice_min << ", w = " << sc.fr_x_slice_w 
        << ", n = " << sc.fr_x_slice_n << "\n";
    return os;
}

//---------------------------------------------------------------------------------------------

std::ostream &operator<<(std::ostream &os, PPSAlignmentConfig c)
{
    os << "* input files\n";
    for (const auto &f : c.inputFiles_)   
        os << "    " << f.c_str() << "\n";
    os << "\n";

    os << "* general info\n";
    os << "    fill = " << c.fill_ << "\n";
    os << "    xangle = " << c.xangle_ << "\n";
    os << "    beta = " << std::fixed << std::setprecision(2) << c.beta_ << "\n";
    os << "    dataset = " << c.dataset_ << "\n\n";

    os << "* dataset already aligned\n";
    os << "    aligned = " << c.aligned_ << "\n\n";

    os << "* alignment parameters\n" << std::setprecision(3);
    for (const auto &p : c.alignmentCorrectionsX_)
        os << "    RP " << p.first << ": de_x = " << p.second << " mm\n";
    os << "\n";

    os << "* cuts\n";
    os << "    n_si = " << c.n_si_ << "\n\n";


    os << "* sector 45\n" << c.sectorConfig45_;
    os << "* sector 56\n" << c.sectorConfig56_ << "\n";

    os << "* matching\n" << std::setprecision(3);
    os << "    reference datasets (" << c.matchingReferenceDatasets_.size() << "):\n";
    for (const auto &ds : c.matchingReferenceDatasets_)
        os << "        " << ds << "\n";

    os << "    shift ranges:\n";
	for (const auto &p : c.matchingShiftRanges_)
        os << "        RP " << p.first << ": sh_min = " << p.second.x_min << ", sh_max = " << p.second.x_max << "\n";

    os << "\n" << "* alignment_x_meth_o\n";
	for (const auto &p : c.alignment_x_meth_o_ranges_)
		os << "    RP " << p.first << ": x_min = " << p.second.x_min << ", x_max = " << p.second.x_max << "\n";

    os << "\n" << "* alignment_x_relative\n";
	for (const auto &p : c.alignment_x_relative_ranges_)
		os << "    RP " << p.first << ": x_min = " << p.second.x_min << ", x_max = " << p.second.x_max << "\n";

    os << "\n" << "* alignment_y\n";
	for (const auto &p : c.alignment_y_ranges_)
		os << "    RP " << p.first << ": x_min = " << p.second.x_min << ", x_max = " << p.second.x_max << "\n";

    return os;
}

TYPELOOKUP_DATA_REG(PPSAlignmentConfig);

#endif 