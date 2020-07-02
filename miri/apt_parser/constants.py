# We store the URI to avoid specifying it everywhere.
ns = {'apt': 'http://www.stsci.edu/JWST/APT',
      "mmrscge": "http://www.stsci.edu/JWST/APT/Template/MiriMRSCrossGratingEngineering",
      "nsmsasd": "http://www.stsci.edu/JWST/APT/Template/NirspecMSAShortDetect",
      "nid": "http://www.stsci.edu/JWST/APT/Template/NirissDark",
      "ncipri": "http://www.stsci.edu/JWST/APT/Template/NircamIprImaging",
      "nif": "http://www.stsci.edu/JWST/APT/Template/NirissFocus",
      "wfscfp": "http://www.stsci.edu/JWST/APT/Template/WfscFinePhasing",
      "nii": "http://www.stsci.edu/JWST/APT/Template/NirissImaging",
      "mmimf": "http://www.stsci.edu/JWST/APT/Template/MiriMimf",
      "mlrs": "http://www.stsci.edu/JWST/APT/Template/MiriLRS",
      "mc": "http://www.stsci.edu/JWST/APT/Template/MiriCoron",
      "md": "http://www.stsci.edu/JWST/APT/Template/MiriDark",
      "nsfgwt": "http://www.stsci.edu/JWST/APT/Template/NirspecFilterGratingWheelTest",
      "wfscga": "http://www.stsci.edu/JWST/APT/Template/WfscGlobalAlignment",
      "nsil": "http://www.stsci.edu/JWST/APT/Template/NirspecInternalLamp",
      "idfu": "http://www.stsci.edu/JWST/APT/Template/IsimDictionaryFileUpdate",
      "mi": "http://www.stsci.edu/JWST/APT/Template/MiriImaging",
      "nsimg": "http://www.stsci.edu/JWST/APT/Template/NirspecImaging",
      "niami": "http://www.stsci.edu/JWST/APT/Template/NirissAmi",
      "niwfss": "http://www.stsci.edu/JWST/APT/Template/NirissWfss",
      "ncif": "http://www.stsci.edu/JWST/APT/Template/NircamInternalFlat",
      "nsfss": "http://www.stsci.edu/JWST/APT/Template/NirspecFixedSlitSpectroscopy",
      "fgsif": "http://www.stsci.edu/JWST/APT/Template/FgsInternalFlat",
      "ncef": "http://www.stsci.edu/JWST/APT/Template/NircamExternalFlat",
      "ncei": "http://www.stsci.edu/JWST/APT/Template/NircamEngineeringImaging",
      "niec": "http://www.stsci.edu/JWST/APT/Template/NirissExternalCalibration",
      "nsmsam": "http://www.stsci.edu/JWST/APT/Template/NirspecMSAMasking",
      "niif": "http://www.stsci.edu/JWST/APT/Template/NirissInternalFlat",
      "ncpili": "http://www.stsci.edu/JWST/APT/Template/NircamPilImaging",
      "nsmsaa": "http://www.stsci.edu/JWST/APT/Template/NirspecMSAAnneal",
      "mcpc": "http://www.stsci.edu/JWST/APT/Template/MiriCpc",
      "fgsec": "http://www.stsci.edu/JWST/APT/Template/FgsExternalCalibration",
      "nsd": "http://www.stsci.edu/JWST/APT/Template/NirspecDark",
      "nsf": "http://www.stsci.edu/JWST/APT/Template/NirspecFocus",
      "ns": "http://www.stsci.edu/JWST/APT/Instrument/Nirspec",
      "nisoss": "http://www.stsci.edu/JWST/APT/Template/NirissSoss",
      "nsmimf": "http://www.stsci.edu/JWST/APT/Template/NirspecMimf",
      "ncts": "http://www.stsci.edu/JWST/APT/Template/NircamTimeSeries",
      "ncd": "http://www.stsci.edu/JWST/APT/Template/NircamDark",
      "mef": "http://www.stsci.edu/JWST/APT/Template/MiriExternalFlat",
      "ncc": "http://www.stsci.edu/JWST/APT/Template/NircamCoron",
      "ncf": "http://www.stsci.edu/JWST/APT/Template/NircamFocus",
      "mmrs": "http://www.stsci.edu/JWST/APT/Template/MiriMRS",
      "nci": "http://www.stsci.edu/JWST/APT/Template/NircamImaging",
      "sk": "http://www.stsci.edu/JWST/APT/Template/StationKeeping",
      "wfscc": "http://www.stsci.edu/JWST/APT/Template/WfscCommissioning",
      "iat": "http://www.stsci.edu/JWST/APT/Template/IsimAsicTuning",
      "sr": "http://www.stsci.edu/JWST/APT/Template/SafeModeRecovery",
      "rtc": "http://www.stsci.edu/JWST/APT/Template/RealtimeCommanding",
      "nsfr": "http://www.stsci.edu/JWST/APT/Template/NirspecFocusReference",
      "ncw": "http://www.stsci.edu/JWST/APT/Template/NircamWheelExercise",
      "nsbots": "http://www.stsci.edu/JWST/APT/Template/NirspecBrightObjectTimeSeries",
      "mann": "http://www.stsci.edu/JWST/APT/Template/MiriAnneal",
      "nsmos": "http://www.stsci.edu/JWST/APT/Template/NirspecMOS",
      "ncgts": "http://www.stsci.edu/JWST/APT/Template/NircamGrismTimeSeries",
      "wfsccp": "http://www.stsci.edu/JWST/APT/Template/WfscCoarsePhasing",
      "nsifus": "http://www.stsci.edu/JWST/APT/Template/NirspecIFUSpectroscopy",
      "msa": "http://www.stsci.edu/JWST/APT/Template/NirspecMSA",
      "ncwfss": "http://www.stsci.edu/JWST/APT/Template/NircamWfss",
      "nifs": "http://www.stsci.edu/JWST/APT/Template/NirissFlatSuite",
      "fgsf": "http://www.stsci.edu/JWST/APT/Template/FgsFocus",
      "po": "http://www.stsci.edu/JWST/APT/Template/PointingOnly"}

# Given the STScI proposal number, return the corresponding MIRI CAR number
MIRI_CAR = {"1024": "MIRI-007", "1027": "MIRI-011", "1028": "MIRI-012",
            "1029": "MIRI-013", "1030": "MIRI-014", "1031": "MIRI-015",
            "1032": "MIRI-016", "1033": "MIRI-017", "1034": "MIRI-018",
            "1037": "MIRI-050", "1038": "MIRI-051", "1039": "MIRI-052",
            "1040": "MIRI-053", "1042": "MIRI-055", "1043": "MIRI-056",
            "1045": "MIRI-058", "1046": "MIRI-060", "1047": "MIRI-061",
            "1048": "MIRI-062", "1049": "MIRI-063", "1050": "MIRI-064",
            "1051": "MIRI-072", "1052": "MIRI-073", "1053": "MIRI-074",
            "1054": "MIRI-075", "1171": "OTE-34", "1172": "OTE-35",
            "1173": "MIRI-001", "1259": "MIRI-076", "1261": "MIRI-077",
            "1406": "MIRI-078", "1011": "MIRI-002", "1012": "MIRI-004",
            "1023": "MIRI-005"}

# Invert previous dict
MIRI_STSCI_ID = {v: k for k, v in MIRI_CAR.items()}

# For some dither patterns, number of points is not necessarily
# specified in the APT file. We then encode for these dither
# pattern how many points there is
DITHER_POINTS = {"4-Point-Sets": 4, "2-Point": 2, "REULEAUX": 12,  # Imager
                 "ALONG_SLIT_NOD": 2, "1-PIXEL_SLIT_SCAN": 45, "2-PIXEL_SLIT_SCAN": 23, "7-PIXEL_SLIT_SCAN": 7,  # LRS
                 "7x3_PIXEL_MAP_CENTER": 21, "7x3_PIXEL_MAP_NOD1": 21, "7x3_PIXEL_MAP_NOD2": 21,  # LRS
                 "SHORT_CROSS_SCAN_NOD_1": 9, "SHORT_CROSS_SCAN_NOD_2": 9, "INTRAPIXEL_SLIT_SCAN_NOD_1": 11,
                 "INTRAPIXEL_SLIT_SCAN_NOD_2": 11,  # LRS
                 "INTRAPIXEL_SLIT_SCAN_CENTER": 11,
                 # LRS. A lot are missing but in the .xsd spec, only slit nod and mapping is referenced
                 "4-POINT-MIRI-F770W-WITH-NIRCam": 4, "4-POINT-MIRI-F1000W-WITH-NIRCam": 4,  # Dual mode with NIRCam
                 "4-POINT-MIRI-F1800W-WITH-NIRCam": 4, "4-POINT-MIRI-F2550W-WITH-NIRCam": 4,  # Dual mode with NIRCam
                 "4-Point": 4,  # DUAL mode Imager/MRS
                 "MRS 2-Point": 2, "MRS 4-Point": 4,  # MRS
                 "5-POINT-SMALL-GRID": 5, "9-POINT-SMALL-GRID": 9,  # Coronograph
                 'CYCLING-MICRO': 16,  # Specific to APT 1028 for PSF super resolution
                 }

SUBARRAY_PIX = {
    "BOTH": 1032 * 1024,  # Artificial subarray for MRS mem and time prediction
    "FULL": 1032 * 1024,
    "BRIGHTSKY": 512 * 512,
    "SUB256": 256 * 256,
    "SUB128": 136 * 128,
    "SUB64": 72 * 64,
    "SLITLESSPRISM": 72 * 416,
    "MASK1065": 288 * 224,
    "MASK1140": 288 * 224,
    "MASK1550": 288 * 224,
    "FOUR_QPM": 288 * 224,
    "MASKLYOT": 320 * 304,
    "LYOT": 320 * 304,
}

DISPERSER = {"SHORT(A)": "SHORT", "MEDIUM(B)": "MEDIUM", "LONG(C)": "LONG"}
