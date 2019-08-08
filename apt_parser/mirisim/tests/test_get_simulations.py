import pytest

sim = pytest.importorskip("apt_parser.mirisim")
import apt_parser.mirisim.constants as c

from apt_parser.utils import assertDictEqual

imsim_testdata = [
    ({'readDetect': 'FAST', 'instrument': 'MIRI', 'subarray': 'BRIGHTSKY', 'filename': '1232.aptx', 'ima_frames': 16,
 'obs_id': '1', 'filter': 'F560W', 'ima_integrations': 18, 'exposures': 1, 'template': 'MiriImaging', 'NDither': 4},
     {"name": "obs1_IMA_F560W", "POP": 'IMA', "ConfigPath": 'IMA_BRIGHTSKY', "Dither": True, "NDither": 4,
      "DitherPat": "ima_recommended_dither.dat", "filter": "F560W", "readDetect": 'BRIGHTSKY',
      "ima_mode": 'FAST', "ima_exposures": 1, "ima_integrations": 18, "ima_frames": 16}
     ),

    ({'readDetect': 'FAST', 'instrument': 'MIRI', 'subarray': 'BRIGHTSKY', 'filename': '1232.aptx', 'ima_frames': 14,
 'obs_id': '1', 'filter': 'F1000W', 'ima_integrations': 16, 'exposures': 1, 'template': 'MiriImaging', 'NDither': 4},
     {"name": "obs1_IMA_F1000W", "POP": 'IMA', "ConfigPath": 'IMA_BRIGHTSKY', "Dither": True, "NDither": 4,
      "DitherPat": "ima_recommended_dither.dat", "filter": "F1000W", "readDetect": 'BRIGHTSKY',
      "ima_mode": 'FAST', "ima_exposures": 1, "ima_integrations": 16, "ima_frames": 14}
     ),

    ({'readDetect': 'FAST', 'instrument': 'MIRI', 'subarray': 'BRIGHTSKY', 'filename': '1232.aptx', 'ima_frames': 14,
 'obs_id': '1', 'filter': 'F1800W', 'ima_integrations': 16, 'exposures': 1, 'template': 'MiriImaging', 'NDither': 4},
     {"name": "obs1_IMA_F1800W", "POP": 'IMA', "ConfigPath": 'IMA_BRIGHTSKY', "Dither": True, "NDither": 4,
      "DitherPat": "ima_recommended_dither.dat", "filter": "F1800W", "readDetect": 'BRIGHTSKY',
      "ima_mode": 'FAST', "ima_exposures": 1, "ima_integrations": 16, "ima_frames": 14}
     ),

    ({'readDetect': 'FAST', 'instrument': 'MIRI', 'subarray': 'BRIGHTSKY', 'filename': '1232.aptx', 'ima_frames': 14,
 'obs_id': '1', 'filter': 'F2550W', 'ima_integrations': 20, 'exposures': 1, 'template': 'MiriImaging', 'NDither': 4},
     {"name": "obs1_IMA_F2550W", "POP": 'IMA', "ConfigPath": 'IMA_BRIGHTSKY', "Dither": True, "NDither": 4,
      "DitherPat": "ima_recommended_dither.dat", "filter": "F2550W", "readDetect": 'BRIGHTSKY',
      "ima_mode": 'FAST', "ima_exposures": 1, "ima_integrations": 20, "ima_frames": 14}
     ),

    ({'instrument': 'MIRI', 'subarray': 'FULL', 'disperser': 'SHORT', 'obs_id': '2', 'filter': 'F560W', 'exposures': 1,
 'LW_frames': 85, 'SW_frames': 85, 'primary_channel': 'ALL', 'filename': '1232.aptx', 'ima_frames': 17,
 'SW_integrations': 3, 'simultaneous_imaging': 'YES', 'ima_integrations': 15, 'readDetect': 'FAST',
 'template': 'MiriMRS', 'NDither': 4, 'detector': 'ALL', 'LW_integrations': 3},
     {"name": "obs2_IMA_F560W", "POP": 'IMA', "ConfigPath": 'IMA_FULL', "Dither": True, "NDither": 4,
      "DitherPat": "mrs_recommended_dither.dat", "filter": "F560W", "readDetect": 'FULL',
      "ima_mode": 'FAST', "ima_exposures": 1, "ima_integrations": 15, "ima_frames": 17}
     ),

    ({'instrument': 'MIRI', 'subarray': 'FULL', 'disperser': 'MEDIUM', 'obs_id': '2', 'filter': 'F770W', 'exposures': 1,
 'LW_frames': 94, 'SW_frames': 94, 'primary_channel': 'ALL', 'filename': '1232.aptx', 'ima_frames': 17,
 'SW_integrations': 3, 'simultaneous_imaging': 'YES', 'ima_integrations': 15, 'readDetect': 'FAST',
 'template': 'MiriMRS', 'NDither': 4, 'detector': 'ALL', 'LW_integrations': 3},
     {"name": "obs2_IMA_F770W", "POP": 'IMA', "ConfigPath": 'IMA_FULL', "Dither": True, "NDither": 4,
      "DitherPat": "mrs_recommended_dither.dat", "filter": "F770W", "readDetect": 'FULL',
      "ima_mode": 'FAST', "ima_exposures": 1, "ima_integrations": 15, "ima_frames": 17}
     ),

    ({'instrument': 'MIRI', 'subarray': 'FULL', 'disperser': 'LONG', 'obs_id': '2', 'filter': 'F1000W', 'exposures': 1,
 'LW_frames': 94, 'SW_frames': 94, 'primary_channel': 'ALL', 'filename': '1232.aptx', 'ima_frames': 17,
 'SW_integrations': 3, 'simultaneous_imaging': 'YES', 'ima_integrations': 15, 'readDetect': 'FAST',
 'template': 'MiriMRS', 'NDither': 4, 'detector': 'ALL', 'LW_integrations': 3},
     {"name": "obs2_IMA_F1000W", "POP": 'IMA', "ConfigPath": 'IMA_FULL', "Dither": True, "NDither": 4,
      "DitherPat": "mrs_recommended_dither.dat", "filter": "F1000W", "readDetect": 'FULL',
      "ima_mode": 'FAST', "ima_exposures": 1, "ima_integrations": 15, "ima_frames": 17}
     ),
]


mrssim_testdata = [
    ({'instrument': 'MIRI', 'subarray': 'FULL', 'disperser': 'SHORT', 'obs_id': '2', 'filter': 'F560W', 'exposures': 1,
 'LW_frames': 85, 'SW_frames': 85, 'primary_channel': 'ALL', 'filename': '1232.aptx', 'ima_frames': 17,
 'SW_integrations': 3, 'simultaneous_imaging': 'YES', 'ima_integrations': 15, 'readDetect': 'FAST',
 'template': 'MiriMRS', 'NDither': 4, 'detector': 'ALL', 'LW_integrations': 3},
     {"name": "obs2_MRS_SHORT", "POP": 'MRS', "ConfigPath": 'MRS_1SHORT', "Dither": True, "NDither": 4,
      "DitherPat": "mrs_recommended_dither.dat", "disperser": "SHORT", "detector": 'BOTH',
      "mrs_mode": 'FAST', "mrs_exposures": 1, "mrs_integrations": 3, "mrs_frames": 85}
     ),

    ({'instrument': 'MIRI', 'subarray': 'FULL', 'disperser': 'MEDIUM', 'obs_id': '2', 'filter': 'F770W', 'exposures': 1,
 'LW_frames': 94, 'SW_frames': 94, 'primary_channel': 'ALL', 'filename': '1232.aptx', 'ima_frames': 17,
 'SW_integrations': 3, 'simultaneous_imaging': 'YES', 'ima_integrations': 15, 'readDetect': 'FAST',
 'template': 'MiriMRS', 'NDither': 4, 'detector': 'ALL', 'LW_integrations': 3},
    {"name": "obs2_MRS_MEDIUM", "POP": 'MRS', "ConfigPath": 'MRS_1MEDIUM', "Dither": True, "NDither": 4,
      "DitherPat": "mrs_recommended_dither.dat", "disperser": "MEDIUM", "detector": 'BOTH',
      "mrs_mode": 'FAST', "mrs_exposures": 1, "mrs_integrations": 3, "mrs_frames": 94}
     ),

    ({'instrument': 'MIRI', 'subarray': 'FULL', 'disperser': 'LONG', 'obs_id': '2', 'filter': 'F1000W', 'exposures': 1,
 'LW_frames': 94, 'SW_frames': 94, 'primary_channel': 'ALL', 'filename': '1232.aptx', 'ima_frames': 17,
 'SW_integrations': 3, 'simultaneous_imaging': 'YES', 'ima_integrations': 15, 'readDetect': 'FAST',
 'template': 'MiriMRS', 'NDither': 4, 'detector': 'ALL', 'LW_integrations': 3},
     {"name": "obs2_MRS_LONG", "POP": 'MRS', "ConfigPath": 'MRS_1LONG', "Dither": True, "NDither": 4,
      "DitherPat": "mrs_recommended_dither.dat", "disperser": "LONG", "detector": 'BOTH',
      "mrs_mode": 'FAST', "mrs_exposures": 1, "mrs_integrations": 3, "mrs_frames": 94}
     ),
]

lrssim_testdata = [
    ({'NDither': 45,  'exposures': 1,  'filename': '1042.aptx',  'ima_frames': 10,  'ima_integrations': 3,
  'instrument': 'MIRI',  'obs_id': '1',  'readDetect': 'FAST',  'subarray': 'FULL',  'template': 'MiriLRS'},
    {"name": "obs1_LRS_SLIT", "POP": 'IMA', "ConfigPath": 'LRS_SLIT', "Dither": True, "NDither": 45,
      "DitherPat": "lrs_recommended_dither.dat", "filter": "P750L", "readDetect": 'FULL',
      "ima_mode": 'FAST', "ima_exposures": 1, "ima_integrations": 3, "ima_frames": 10}
     ),
    ( {'NDither': 1,  'exposures': 1,  'filename': '1042.aptx',  'ima_frames': 125,  'ima_integrations': 4,
  'instrument': 'MIRI',  'obs_id': '5',  'readDetect': 'FAST',  'subarray': 'SLITLESSPRISM',  'template': 'MiriLRS'},
    {"name": "obs5_LRS_SLITLESS", "POP": 'IMA', "ConfigPath": 'LRS_SLITLESS', "Dither": False, "NDither": 1,
      "DitherPat": "lrs_recommended_dither.dat", "filter": "P750L", "readDetect": 'SLITLESSPRISM',
      "ima_mode": 'FAST', "ima_exposures": 1, "ima_integrations": 4, "ima_frames": 125}
     ),
]





@pytest.mark.parametrize("observation,expected_simulation", imsim_testdata)
def test_get_imsim(observation, expected_simulation):

    simulation = sim.get_imsim_parameters(observation)

    # Only mandatory values are stored for the test. The others come from the default dictionnary, dynamically
    ref_sim = c.MIRISIM_DEFAULT.copy()
    ref_sim.update(expected_simulation)  # Change only the paramters given in the simulation, all others are the default one.

    assertDictEqual(simulation, ref_sim)


@pytest.mark.parametrize("observation,expected_simulation", lrssim_testdata)
def test_get_lrssim(observation, expected_simulation):

    simulation = sim.get_lrs_parameters(observation)

    # Only mandatory values are stored for the test. The others come from the default dictionnary, dynamically
    ref_sim = c.MIRISIM_DEFAULT.copy()
    ref_sim.update(expected_simulation)  # Change only the paramters given in the simulation, all others are the default one.

    assertDictEqual(simulation, ref_sim)


@pytest.mark.parametrize("observation,expected_simulation", mrssim_testdata)
def test_get_maisie(observation, expected_simulation):

    simulation = sim.get_mrs_parameters(observation)

    # Only mandatory values are stored for the test. The others come from the default dictionnary, dynamically
    ref_sim = c.MIRISIM_DEFAULT.copy()
    ref_sim.update(expected_simulation)  # Change only the paramters given in the simulation, all others are the default one.

    assertDictEqual(simulation, ref_sim)
