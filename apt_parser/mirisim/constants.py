MIRISIM_DEFAULT = {
"name":"Default Simulation", "rel_obsdate":0.0, "scene":"scene.ini",
"POP":'IMA', "ConfigPath":'IMA_FULL', "Dither":False, "StartInd":1, "NDither":2,
"DitherPat":"ima_recommended_dither.dat",
"filter":"F560W", "readDetect":'FULL', "ima_mode":'FAST', "ima_exposures":1, "ima_integrations":1, "ima_frames":1,
"disperser":'SHORT', "detector":'SW', "mrs_mode":'FAST', "mrs_exposures":1, "mrs_integrations":1, "mrs_frames":1
}

PRIMARY_CHANNEL = {"ALL":1, "CHANNEL1":1, "CHANNEL2":2, "CHANNEL3":3, "CHANNEL4":4}

NON_SUPPORTED_TEMPLATES = ["MiriDark", "MiriExternalFlat", "MiriCoron", "MiriCpc", 'MiriMRSCrossGratingEngineering']