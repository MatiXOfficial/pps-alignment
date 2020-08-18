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

#include <iostream>
#include <cmath>
#include <iomanip>

// -------------------------------- PPSAlignmentConfig getters --------------------------------

std::vector<std::string> PPSAlignmentConfig::sequence() const { return sequence_; }

SectorConfig PPSAlignmentConfig::sectorConfig45() const { return sectorConfig45_; }
SectorConfig PPSAlignmentConfig::sectorConfig56() const { return sectorConfig56_; }

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

std::vector<std::string> PPSAlignmentConfig::matchingReferenceDatasets() const 
{ 
	return matchingReferenceDatasets_; 
}
std::map<unsigned int, SelectionRange> PPSAlignmentConfig::matchingShiftRanges() const 
{ 
	return matchingShiftRanges_; 
}

std::map<unsigned int, double> PPSAlignmentConfig::yMaxFit() const 
{ 
	return yMaxFit_; 
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

void PPSAlignmentConfig::setSequence(std::vector<std::string> &sequence) { sequence_ = sequence; }

void PPSAlignmentConfig::setSectorConfig45(SectorConfig &sectorConfig45) { sectorConfig45_ = sectorConfig45; }
void PPSAlignmentConfig::setSectorConfig56(SectorConfig &sectorConfig56) { sectorConfig56_ = sectorConfig56; }

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

void PPSAlignmentConfig::setMatchingReferenceDatasets(std::vector<std::string> &matchingReferenceDatasets)
{
	matchingReferenceDatasets_ = matchingReferenceDatasets;
}
void PPSAlignmentConfig::setMatchingShiftRanges(std::map<unsigned int, SelectionRange> &matchingShiftRanges)
{
	matchingShiftRanges_ = matchingShiftRanges;
}

void PPSAlignmentConfig::setYMaxFit(std::map<unsigned int, double> &yMaxFit)
{
	yMaxFit_ = yMaxFit;
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

std::ostream &operator<<(std::ostream &os, RPConfig &rc)
{
	os << std::fixed << std::setprecision(3);
	os << "    " << rc.name << ", id = " << rc.id << ", position = " << rc.position << ":\n";
	os << "        slope = " << rc.slope << ", sh_x = " << rc.sh_x;

	return os;
}

std::ostream &operator<<(std::ostream &os, SectorConfig &sc)
{
	// to be adjusted
	os << std::fixed << std::setprecision(3);
	os << sc.name << ":\n";
	os << sc.rp_N << "\n" << sc.rp_F << "\n";
	os << "    slope = " << sc.slope << "\n";
	os << "    cut_h: apply = " << sc.cut_h_apply << ", a = " << sc.cut_h_a << ", c = " 
		<< sc.cut_h_c << ", si = " << sc.cut_h_si << "\n";
	os << "    cut_v: apply = " << sc.cut_v_apply << ", a = " << sc.cut_v_a << ", c = " 
		<< sc.cut_v_c << ", si = " << sc.cut_v_si << "\n";
	os << std::setprecision(2);
	os << "    x slices, nr: min = " << sc.nr_x_slice_min << ", w = " << sc.nr_x_slice_w 
		<< ", n = " << sc.nr_x_slice_n << "\n";
	os << "    x slices, fr: min = " << sc.fr_x_slice_min << ", w = " << sc.fr_x_slice_w 
		<< ", n = " << sc.fr_x_slice_n;

	return os;
}

std::ostream &operator<<(std::ostream &os, PPSAlignmentConfig c)
{
	os << "* sequence\n";
	for (unsigned int i = 0; i < c.sequence_.size(); i++)
	{
		os << "    " << i + 1 << ": " << c.sequence_[i] << "\n";
	}
	os << "\n";

	os << "* " << c.sectorConfig45_ << "\n\n";
	os << "* " << c.sectorConfig56_ << "\n\n";

	os << "* alignment corrections\n" << std::setprecision(3);
	for (const auto &p : c.alignmentCorrectionsX_)
	{
		os << "    RP " << p.first << ": de_x = " << p.second << " mm, de_y = " 
		<< c.alignmentCorrectionsY_[p.first] << " mm\n";
	}
	os << "\n";

	os << "* dataset already aligned\n";
	os << "    aligned = " << c.aligned_ << "\n\n";

	os << "* cuts\n";
	os << "    n_si = " << c.n_si_ << "\n\n";

	os << "* matching\n" << std::setprecision(3);
	os << "    reference datasets (" << c.matchingReferenceDatasets_.size() << "):\n";
	for (const auto &ds : c.matchingReferenceDatasets_)
		os << "        " << ds << "\n";

	os << "    shift ranges:\n";
	for (const auto &p : c.matchingShiftRanges_)
		os << "        RP " << p.first << ": sh_min = " << p.second.x_min << ", sh_max = " << p.second.x_max << "\n";
		
	os << "\n" << "* y_max_fit\n";
	for (const auto &p : c.yMaxFit_)
		os << "    RP " << p.first << ": " << p.second << "\n";

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