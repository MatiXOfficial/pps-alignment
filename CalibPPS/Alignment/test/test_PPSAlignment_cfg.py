import FWCore.ParameterSet.Config as cms

process = cms.Process('testAlignment')

process.load('FWCore.MessageService.MessageLogger_cfi')

process.MessageLogger = cms.Service("MessageLogger",
	destinations = cms.untracked.vstring('alignment_out', 
	                                     'alignment_log', 
	                                     'cout'
										 ),
	categories = cms.untracked.vstring('x_alignment_results',
                                       'x_alignment_relative_results',
                                       'y_alignment_results', 
									   ),
	alignment_out = cms.untracked.PSet(
		threshold = cms.untracked.string("INFO"),
		INFO = cms.untracked.PSet(limit = cms.untracked.int32(0)),
		x_alignment_results = cms.untracked.PSet(limit = cms.untracked.int32(1000)),
		x_alignment_relative_results = cms.untracked.PSet(limit = cms.untracked.int32(1000)),
		y_alignment_results = cms.untracked.PSet(limit = cms.untracked.int32(1000))
	),
	alignment_log = cms.untracked.PSet(
		threshold = cms.untracked.string("INFO")
	),
	cout = cms.untracked.PSet(
		threshold = cms.untracked.string('WARNING')
	)
)

# load DQM framework
process.load("DQM.Integration.config.environment_cfi")
process.dqmEnv.subSystemFolder = "CalibPPS"
process.dqmEnv.eventInfoFolder = "EventInfo"
process.dqmSaver.path = ""
process.dqmSaver.tag = "CalibPPS"

process.source = cms.Source("PoolSource",
	fileNames = cms.untracked.vstring('root://eostotem.cern.ch//eos/cms/store/group/phys_pps/reconstruction/2018/physics_runs/rec-hit-version1/fill7334_xangle160_beta0.30_EGamma.root')
)

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

process.ppsAlignmentConfigESSource = cms.ESSource("PPSAlignmentConfigESSource",
	sequence = cms.vstring(
		"x alignment",
		"x alignment relative",
		"y alignment"
	),
	
	sector_45 = cms.PSet(
		rp_N = cms.PSet(
			slope = cms.double(0.19),
			sh_x = cms.double(-3.6)
		),
		rp_F = cms.PSet(
			slope = cms.double(0.19),
			sh_x = cms.double(-42.)
		),
		slope = cms.double(0.006),

		cut_h_apply = cms.bool(True),
		cut_h_a = cms.double(-1),
		cut_h_c = cms.double(-38.55 + 0.57 - 0.08),
		cut_h_si = cms.double(0.2),

		cut_v_apply = cms.bool(True),
		cut_v_a = cms.double(-1.07),
		cut_v_c = cms.double(1.63 - 2.15 + 0.25),
		cut_v_si = cms.double(0.15),

		nr_x_slice_min = cms.double(7),
		nr_x_slice_max = cms.double(19),
		nr_x_slice_w = cms.double(0.2),

		fr_x_slice_min = cms.double(46),
		fr_x_slice_max = cms.double(58),
		fr_x_slice_w = cms.double(0.2),
	),

	sector_56 = cms.PSet(
		rp_N = cms.PSet(
			slope = cms.double(0.40),
			sh_x = cms.double(-2.8)
		),
		rp_F = cms.PSet(
			slope = cms.double(0.39),
			sh_x = cms.double(-41.9)
		),
		slope = cms.double(-0.015),

		cut_h_apply = cms.bool(True),
		cut_h_a = cms.double(-1),
		cut_h_c = cms.double(-39.26 + 0.33),
		cut_h_si = cms.double(0.2),

		cut_v_apply = cms.bool(True),
		cut_v_a = cms.double(-1.07),
		cut_v_c = cms.double(1.49 - 1.80),
		cut_v_si = cms.double(0.15),

		nr_x_slice_min = cms.double(6),
		nr_x_slice_max = cms.double(17.),
		nr_x_slice_w = cms.double(0.2),

		fr_x_slice_min = cms.double(45),
		fr_x_slice_max = cms.double(57.),
		fr_x_slice_w = cms.double(0.2),
	),

	aligned = cms.bool(False),
	n_si = cms.double(4.),

	matching = cms.PSet(
		reference_datasets = cms.vstring("/eos/user/k/kocotm/Documents/Alignment_data/alig-version3/fill_6554/xangle_130_beta_0.30/DS1/distributions.root"),

		rp_L_2_F = cms.PSet(
			sh_min = cms.double(-43),
			sh_max = cms.double(-41)
		),
		rp_L_1_F = cms.PSet(
			sh_min = cms.double(-4.2),
			sh_max = cms.double(-2.4)
		),
		rp_R_1_F = cms.PSet(
			sh_min = cms.double(-3.6),
			sh_max = cms.double(-1.8)
		),
		rp_R_2_F = cms.PSet(
			sh_min = cms.double(-43.2),
			sh_max = cms.double(-41.2)
		)
	),
	
	y_max_fit = cms.PSet(
		rp_L_2_F = cms.double(7.5),
		rp_L_1_F = cms.double(7.8),
		rp_R_1_F = cms.double(7.4),
		rp_R_2_F = cms.double(8.0)
	),

	x_alignment_meth_o = cms.PSet(
	rp_L_2_F = cms.PSet(
		x_min = cms.double(47.),
		x_max = cms.double(56.5),
	),
	rp_L_1_F = cms.PSet(
		x_min = cms.double(9.),
		x_max = cms.double(18.5),
	),
	rp_R_1_F = cms.PSet(
		x_min = cms.double(7.),
		x_max = cms.double(15.),
	),
	rp_R_2_F = cms.PSet(
		x_min = cms.double(46.),
		x_max = cms.double(54.),
	)
	),

	x_alignment_relative = cms.PSet(
	rp_L_2_F = cms.PSet(
		x_min = cms.double(0.),
		x_max = cms.double(0.),
	),
	rp_L_1_F = cms.PSet(
		x_min = cms.double(7.5),
		x_max = cms.double(12.),
	),
	rp_R_1_F = cms.PSet(
		x_min = cms.double(6.),
		x_max = cms.double(10.),
	),
	rp_R_2_F = cms.PSet(
		x_min = cms.double(0.),
		x_max = cms.double(0.),
	)
	),

	y_alignment = cms.PSet(
	rp_L_2_F = cms.PSet(
		x_min = cms.double(44.5),
		x_max = cms.double(49.0),
	),
	rp_L_1_F = cms.PSet(
		x_min = cms.double(6.7),
		x_max = cms.double(11.0),
	),
	rp_R_1_F = cms.PSet(
		x_min = cms.double(5.9),
		x_max = cms.double(10.0),
	),
	rp_R_2_F = cms.PSet(
		x_min = cms.double(44.5),
		x_max = cms.double(49.0),
	)
	)
)

process.ppsAlignmentConfigESSource_reference = cms.ESSource("PPSAlignmentConfigESSource",
	label = cms.string('reference'),

	aligned = cms.bool(True),
	
	y_max_fit = cms.PSet(
		rp_L_2_F = cms.double(3.5),
		rp_L_1_F = cms.double(4.5),
		rp_R_1_F = cms.double(5.5),
		rp_R_2_F = cms.double(4.8)
	),

	x_alignment_meth_o = cms.PSet(
	rp_L_2_F = cms.PSet(
		x_min = cms.double(5.),
		x_max = cms.double(15.),
	),
	rp_L_1_F = cms.PSet(
		x_min = cms.double(5.),
		x_max = cms.double(15.),
	),
	rp_R_1_F = cms.PSet(
		x_min = cms.double(4.),
		x_max = cms.double(12.),
	),
	rp_R_2_F = cms.PSet(
		x_min = cms.double(4.),
		x_max = cms.double(12.),
	)
	),

	binning = cms.PSet(
		pixel_x_offset = cms.double(0.)
	)
)

process.load("CalibPPS.Alignment.ppsAlignmentWorker_cfi")
process.load("CalibPPS.Alignment.ppsAlignmentHarvester_cfi")

process.path = cms.Path(
  	process.ppsAlignmentWorker
  	* process.ppsAlignmentHarvester
)

process.end_path = cms.EndPath(
  process.dqmEnv +
  process.dqmSaver
)

process.schedule = cms.Schedule(
  process.path,
  process.end_path
)