/****************************************************************************
 *
 *  CondFormats/PPSObjects/interface/PPSAlignmentConfig.cc
 *
 *  Description : Alignment parameters
 *
 * Authors:
 *  - Jan Kašpar
 *  - Mateusz Kocot
 *
 ****************************************************************************/

#include "FWCore/Utilities/interface/typelookup.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "CondFormats/PPSObjects/interface/PPSAlignmentConfig.h"
TYPELOOKUP_DATA_REG(PPSAlignmentConfig);

#include <cmath>
#include <iomanip>

// -------------------------------- SelectionRange --------------------------------

SelectionRange::SelectionRange(double x_min, double x_max) : x_min(x_min), x_max(x_max)
{}

// -------------------------------- PPSAlignmentConfig getters --------------------------------

unsigned int PPSAlignmentConfig::fill() const { return fill_; }
unsigned int PPSAlignmentConfig::xangle() const { return xangle_; }
double PPSAlignmentConfig::beta() const { return beta_; }
std::string PPSAlignmentConfig::dataset() const { return dataset_; }

std::map<unsigned int, std::string> PPSAlignmentConfig::rpTags() const { return rpTags_; }

std::vector<std::string> PPSAlignmentConfig::inputFiles() const { return inputFiles_; }

std::map<unsigned int, double> PPSAlignmentConfig::alignmentCorrectionsX() const 
{ 
    return alignmentCorrectionsX_; 
}
std::map<unsigned int, double> PPSAlignmentConfig::alignmentCorrectionsY() const 
{ 
    return alignmentCorrectionsY_;
}

bool PPSAlignmentConfig::aligned() const { return aligned_; }

double PPSAlignmentConfig::n_si() const { return n_si_; }

SectorConfig PPSAlignmentConfig::sectorConfig45() const { return sectorConfig45_; }
SectorConfig PPSAlignmentConfig::sectorConfig56() const { return sectorConfig56_; }

std::vector<std::string> PPSAlignmentConfig::matchingReferenceDatasets() const 
{ 
    return matchingReferenceDatasets_; 
}
std::map<unsigned int, SelectionRange> PPSAlignmentConfig::matchingShiftRanges() const 
{ 
    return matchingShiftRanges_; 
}

std::map<unsigned int, SelectionRange> PPSAlignmentConfig::alignment_x_meth_o_ranges() const 
{ 
    return alignment_x_meth_o_ranges_; 
}
std::map<unsigned int, SelectionRange> PPSAlignmentConfig::alignment_x_relative_ranges() const 
{ 
    return alignment_x_relative_ranges_; 
}

std::map<unsigned int, SelectionRange> PPSAlignmentConfig::alignment_y_ranges() const 
{ 
    return alignment_y_ranges_; 
}

// -------------------------------- PPSAlignmentConfig setters --------------------------------

void PPSAlignmentConfig::setFill(unsigned int fill) { fill_ = fill; }
void PPSAlignmentConfig::setXangle(unsigned int xangle) { xangle_ = xangle; }
void PPSAlignmentConfig::setBeta(double beta) { beta_ = beta; }
void PPSAlignmentConfig::setDataset(std::string &dataset) { dataset_ = dataset; }

void PPSAlignmentConfig::setRpTags(std::map<unsigned int, std::string> &rpTags) { rpTags_ = rpTags; }

void PPSAlignmentConfig::setInputFiles(std::vector<std::string> &inputFiles) { inputFiles_ = inputFiles; }

void PPSAlignmentConfig::setAlignmentCorrectionsX(std::map<unsigned int, double> &alignmentCorrectionsX) 
{ 
    alignmentCorrectionsX_ = alignmentCorrectionsX; 
}
void PPSAlignmentConfig::setAlignmentCorrectionsY(std::map<unsigned int, double> &alignmentCorrectionsY) 
{ 
    alignmentCorrectionsY_ = alignmentCorrectionsY; 
}

void PPSAlignmentConfig::setAligned(bool aligned) { aligned_ = aligned; }

void PPSAlignmentConfig::setN_si(double n_si) { n_si_ = n_si; }

void PPSAlignmentConfig::setSectorConfig45(SectorConfig &sectorConfig45) { sectorConfig45_ = sectorConfig45; }
void PPSAlignmentConfig::setSectorConfig56(SectorConfig &sectorConfig56) { sectorConfig56_ = sectorConfig56; }

void PPSAlignmentConfig::setMatchingReferenceDatasets(std::vector<std::string> &matchingReferenceDatasets)
{
    matchingReferenceDatasets_ = matchingReferenceDatasets;
}
void PPSAlignmentConfig::setMatchingShiftRanges(std::map<unsigned int, SelectionRange> &matchingShiftRanges)
{
    matchingShiftRanges_ = matchingShiftRanges;
}

void PPSAlignmentConfig::setAlignment_x_meth_o_ranges(std::map<unsigned int, SelectionRange> &alignment_x_meth_o_ranges)
{
    alignment_x_meth_o_ranges_ = alignment_x_meth_o_ranges;
}
void PPSAlignmentConfig::setAlignment_x_relative_ranges(std::map<unsigned int, SelectionRange> &alignment_x_relative_ranges)
{
    alignment_x_relative_ranges_ = alignment_x_relative_ranges;   
}

void PPSAlignmentConfig::setAlignment_y_ranges(std::map<unsigned int, SelectionRange> &alignment_y_ranges)
{
    alignment_y_ranges_ = alignment_y_ranges;
}

// -------------------------------- << operators --------------------------------

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