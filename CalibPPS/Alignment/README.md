# PPSAlignmentWorker
## Config example:
`CalibPPS/Alignment/python/ppsAlignmentWorker_cfi.py`
## Parameters:
| Name        | Type           | Description                                                              |
|-------------|----------------|--------------------------------------------------------------------------|
| `tagTracks` | `cms.InputTag` |                                                                          |
| `folder`    | `cms.string`   | Should be the same as the `folder` parameter in DQM configuration.       |
| `label`     | `cms.string`   | label for EventSetup                                                     |
| `debug`     | `cms.bool`     | When set to `True`, the worker will produce some extra debug histograms. |

# PPSAlignmentHarvester
## Config example:
`CalibPPS/Alignment/python/ppsAlignmentHarvester_cfi.py`
## Parameters:
| Name     | Type         | Description                                                                                     |
|----------|--------------|-------------------------------------------------------------------------------------------------|
| `folder` | `cms.string` | Should be the same as the `folder` parameter in DQM configuration.                              |
| `debug`  | `cms.bool`   | When set to `True`, the harvester harvester will produce an extra ROOT   file with debug plots. |

# Event Setup
Default values come from the `fillDescriptions` method in `CalibPPS/ESProducers/plugins/PPSAlignmentConfigESSource.cc`.
## Main parameters:
| Name                   | Type          | Default       | Description                                                                                                  |
|------------------------|---------------|---------------|--------------------------------------------------------------------------------------------------------------|
| `debug`                | `cms.bool`    | `False`       | When set to `True`, the ESProducer will produce an extra ROOT file with   debug plots (from reference run).  |
| `label`                | `cms.string`  | `""`          | label to distinguish reference and test fill configs. Should be set   either to `""` (test) or `"reference"` |
| `sequence`             | `cms.vstring` | empty vector  | Determines order of the alignment methods: `"x_alignemnt"`,   `"x_alignment_relative"`, `"y_alignment"`.     |
| `sector_45`            | `cms.PSet`    | details below | Configuration of sector 45. Details below                                                                    |
| `sector_56`            | `cms.PSet`    | details below | Configuration of sector 56. Details below                                                                    |
| `x_ali_sh_step`        | `cms.double`  | `0.01`        | Step for x alignment algorithm                                                                               |
| `y_mode_sys_unc`       | `cms.double`  | `0.03`        | Squared is an element of y mode uncertainity in y alignment (harvester -   line 623).                        |
| `chiSqThreshold`       | `cms.double`  | `50.`         | Chi-square threshold of y mode (harvester - line 626)                                                        |
| `y_mode_unc_max_valid` | `cms.double`  | `5.`          | Maximal valid y mode uncertainity (harvester - line 628)                                                     |
| `y_mode_max_valid`     | `cms.double`  | `20.`         | Maximal valid y mode (harvester - line 629)                                                                  |
| `max_RP_tracks_size`   | `cms.uint32`  | `2.`          | Maximal tracksUp or tracksDw size to avoid crowded events (worker - line   243)                              |
| `n_si`                 | `cms.double`  | `4.`          | Element of determing checking whether the cuts passed (worker - line 268,   272)                             |
| `matching`             | `cms.PSet`    | details below | Reference dataset parameters. Details below                                                                  |
| `x_alignment_meth_o`   | `cms.PSet`    | details below | X alignment parameters. Details below                                                                        |
| `x_alignment_relative` | `cms.PSet`    | details below | Relative x alignment parameters. Details below                                                               |
| `y_alignment`          | `cms.PSet`    | details below | Y alignment parameters. Details below                                                                        |
| `binning`              | `cms.PSet`    | details below | Binning parameters. Details below                                                                            |