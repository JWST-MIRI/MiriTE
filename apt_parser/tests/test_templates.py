from miri.apt_parser import templates
from miri.apt_parser.utils import assertListDictEqual, read_fake_xml


def test_MiriCoron():
    """
    Test of MiriCoron template

    Sample from APT STScI ID 1406 observation 1
    """

    sample_xml = """    <mc:MiriCoron>
        <mc:AcqTargetID>Same Target as Observation</mc:AcqTargetID>
        <mc:AcqFilter>FND</mc:AcqFilter>
        <mc:AcqReadoutPattern>FAST</mc:AcqReadoutPattern>
        <mc:AcqGroups>275</mc:AcqGroups>
        <mc:AcqEtcId></mc:AcqEtcId>
        <mc:AcqQuadrant>1</mc:AcqQuadrant>
        <mc:SecondExposure>YES</mc:SecondExposure>
        <mc:Dither>5-POINT-SMALL-GRID</mc:Dither>
        <mc:Filters>
            <mc:FilterConfig>
                <mc:Mask>FOUR_QPM</mc:Mask>
                <mc:Filter>F1065C</mc:Filter>
                <mc:Exposures>1</mc:Exposures>
                <mc:ReadoutPattern>FAST</mc:ReadoutPattern>
                <mc:Groups>80</mc:Groups>
                <mc:Integrations>10</mc:Integrations>
                <mc:EtcId></mc:EtcId>
            </mc:FilterConfig>
        </mc:Filters>
        <mc:PsfReference>true</mc:PsfReference>
        <mc:PsfJustification>false</mc:PsfJustification>
    </mc:MiriCoron>"""

    tree = read_fake_xml(sample_xml)

    metadata = {'obs_id': '1', 'filename': 'apt/1406.aptx', 'instrument': 'MIRI'}

    obj = templates.MiriCoron(tree, metadata=metadata)
    output_dict = obj.getobs()

    expected_dict = [{'readDetect': 'FAST', 'filter': 'F1065C', 'NDither': 5, 'subarray': 'FOUR_QPM',
                     'obs_id': '1', 'ima_frames': 80, 'filename': 'apt/1406.aptx', 'exposures': 1,
                     'template': 'MiriCoron', 'ima_integrations': 10, 'instrument': 'MIRI'}]

    assertListDictEqual(output_dict, expected_dict)


def test_MiriCpc():
    """
    Test of MiriCpc template

    Sample from APT STScI ID 1045 observation 1
    """

    sample_xml = """    <mcpc:MiriCpc>
        <mcpc:Subarray>MASK1065</mcpc:Subarray>
        <mcpc:Filters>
            <mcpc:FilterConfig>
                <mcpc:Exposures>1</mcpc:Exposures>
                <mcpc:Filter>F1065C</mcpc:Filter>
                <mcpc:ReadoutPattern>FAST</mcpc:ReadoutPattern>
                <mcpc:Groups>100</mcpc:Groups>
                <mcpc:Integrations>1</mcpc:Integrations>
                <mcpc:EtcId></mcpc:EtcId>
            </mcpc:FilterConfig>
        </mcpc:Filters>
    </mcpc:MiriCpc>"""

    tree = read_fake_xml(sample_xml)

    metadata = {'obs_id': '1', 'filename': 'apt/1045.aptx', 'instrument': 'MIRI'}

    obj = templates.MiriCpc(tree, metadata=metadata)
    output_dict = obj.getobs()

    expected_dict = [{'filter': 'F1065C', 'ima_frames': 100, 'obs_id': '1', 'filename': 'apt/1045.aptx',
                     'ima_integrations': 1, 'instrument': 'MIRI', 'exposures': 1, 'readDetect': 'FAST',
                     'template': 'MiriCpc', 'NDither': 1, 'subarray': 'MASK1065', 'instrument': 'MIRI'}]

    assertListDictEqual(output_dict, expected_dict)


def test_MiriDark():
    """
    Test of MiriDark template

    Sample from APT STScI ID 1046 observation 1
    """

    sample_xml = """    <md:MiriDark>
        <md:Detector>IMAGER</md:Detector>
        <md:TestPattern>NONE</md:TestPattern>
        <md:Subarray>FULL</md:Subarray>
        <md:Filters>
            <md:FilterConfig>
                <md:ReadoutPattern>SLOW</md:ReadoutPattern>
                <md:Groups>6</md:Groups>
                <md:Integrations>3</md:Integrations>
                <md:EtcId></md:EtcId>
                <md:Exposures>2</md:Exposures>
            </md:FilterConfig>
        </md:Filters>
    </md:MiriDark>"""

    tree = read_fake_xml(sample_xml)

    metadata = {'obs_id': '1', 'filename': 'apt/1046.aptx', 'instrument': 'MIRI'}

    obj = templates.MiriDark(tree, metadata=metadata)
    output_dict = obj.getobs()

    expected_dict = [{'ima_frames': 6, 'ima_integrations': 3, 'obs_id': '1', 'filename': 'apt/1046.aptx',
                      'exposures': 2, 'subarray': 'FULL', 'instrument': 'MIRI', 'NDither': 1,
                        'readDetect': 'SLOW', 'template': 'MiriDark', 'detector': 'IMAGER'}]

    assertListDictEqual(output_dict, expected_dict)


def test_MiriExternalFlat():
    """
    Test of MiriExternalFlat template

    Sample from APT STScI ID 1027 observation 9
    """

    sample_xml = """    <mef:MiriExternalFlat>
        <mef:Detector>IMAGER</mef:Detector>
        <mef:Dither>true</mef:Dither>
        <mef:DitherSpec>
            <mef:DitherType>4-Point-Sets</mef:DitherType>
            <mef:StartingSet>1</mef:StartingSet>
            <mef:NumberOfSets>4</mef:NumberOfSets>
            <mef:Direction>POSITIVE</mef:Direction>
            <mef:OptimizeFor>POINT SOURCE</mef:OptimizeFor>
            <mef:PatternSize>DEFAULT</mef:PatternSize>
        </mef:DitherSpec>
        <mef:Lamp>OFF THEN ON</mef:Lamp>
        <mef:Subarray>FULL</mef:Subarray>
        <mef:ExposureList>
            <mef:Exposure>
                <mef:Exposures>1</mef:Exposures>
                <mef:Filter>F2550WR</mef:Filter>
                <mef:ReadoutPattern>FAST</mef:ReadoutPattern>
                <mef:Groups>4</mef:Groups>
                <mef:Integrations>1</mef:Integrations>
                <mef:EtcId></mef:EtcId>
                <mef:GroupsShadow>15</mef:GroupsShadow>
                <mef:IntegrationsShadow>1</mef:IntegrationsShadow>
                <mef:EtcIdShadow></mef:EtcIdShadow>
            </mef:Exposure>
        </mef:ExposureList>
    </mef:MiriExternalFlat>"""

    tree = read_fake_xml(sample_xml)

    metadata = {'obs_id': '9', 'filename': 'apt/1027.aptx', 'instrument': 'MIRI'}

    obj = templates.MiriExternalFlat(tree, metadata=metadata)
    output_dict = obj.getobs()

    expected_dict = [{'ima_frames': 4, 'ima_integrations': 1, 'obs_id': '9', 'filename': 'apt/1027.aptx',
                     'dither_type': 'true', 'exposures': 1, 'subarray': 'FULL', 'filter': 'F2550WR', 'instrument': 'MIRI', 'NDither': 16,
                     'readDetect': 'FAST', 'template': 'MiriExternalFlat', 'detector': 'IMAGER'}]

    assertListDictEqual(output_dict, expected_dict)


def test_MiriImaging():
    """
    Test of MiriImaging template

    Sample from APT STScI ID 1052 observation 1
    """

    sample_xml = """    <mi:MiriImaging>
        <mi:Subarray>FULL</mi:Subarray>
        <mi:Dithers>
            <mi:DitherSpecification>
                <mi:DitherType>CYCLING</mi:DitherType>
                <mi:StartingPoint>1</mi:StartingPoint>
                <mi:NumberOfPoints>3</mi:NumberOfPoints>
                <mi:StartingSet>1</mi:StartingSet>
                <mi:NumberOfSets>1</mi:NumberOfSets>
                <mi:PatternSize>DEFAULT</mi:PatternSize>
            </mi:DitherSpecification>
        </mi:Dithers>
        <mi:Filters>
            <mi:FilterConfig>
                <mi:Exposures>1</mi:Exposures>
                <mi:Filter>F770W</mi:Filter>
                <mi:ReadoutPattern>FAST</mi:ReadoutPattern>
                <mi:Groups>100</mi:Groups>
                <mi:Integrations>1</mi:Integrations>
                <mi:EtcId></mi:EtcId>
                <mi:Dither>Dither 1</mi:Dither>
            </mi:FilterConfig>
            <mi:FilterConfig>
                <mi:Exposures>1</mi:Exposures>
                <mi:Filter>F1280W</mi:Filter>
                <mi:ReadoutPattern>FAST</mi:ReadoutPattern>
                <mi:Groups>100</mi:Groups>
                <mi:Integrations>1</mi:Integrations>
                <mi:EtcId></mi:EtcId>
                <mi:Dither>Dither 1</mi:Dither>
            </mi:FilterConfig>
            <mi:FilterConfig>
                <mi:Exposures>1</mi:Exposures>
                <mi:Filter>F1500W</mi:Filter>
                <mi:ReadoutPattern>FAST</mi:ReadoutPattern>
                <mi:Groups>100</mi:Groups>
                <mi:Integrations>1</mi:Integrations>
                <mi:EtcId></mi:EtcId>
                <mi:Dither>Dither 1</mi:Dither>
            </mi:FilterConfig>
        </mi:Filters>
    </mi:MiriImaging>"""

    tree = read_fake_xml(sample_xml)

    metadata = {'obs_id': '1', 'filename': 'apt/1052.aptx', 'instrument': 'MIRI'}

    obj = templates.MiriImaging(tree, metadata=metadata)
    output_dict = obj.getobs()

    expected_dict = [{'filter': 'F770W', 'ima_frames': 100, 'ima_integrations': 1, 'obs_id': '1', 'filename': 'apt/1052.aptx',
                     'exposures': 1, 'readDetect': 'FAST', 'instrument': 'MIRI', 'template': 'MiriImaging',
                     'NDither': 3, 'subarray': 'FULL'},
                     {'filter': 'F1280W', 'ima_frames': 100, 'ima_integrations': 1, 'obs_id': '1', 'filename': 'apt/1052.aptx',
                      'exposures': 1, 'readDetect': 'FAST', 'instrument': 'MIRI', 'template': 'MiriImaging',
                      'NDither': 3, 'subarray': 'FULL'},
                     {'filter': 'F1500W', 'ima_frames': 100, 'ima_integrations': 1, 'obs_id': '1', 'filename': 'apt/1052.aptx',
                      'exposures': 1, 'readDetect': 'FAST', 'instrument': 'MIRI', 'template': 'MiriImaging',
                      'NDither': 3, 'subarray': 'FULL'}
                     ]

    assertListDictEqual(output_dict, expected_dict)


def test_MiriLRS():
    """
    Test of MiriLRS template

    Sample from APT STScI ID 1042 observation 1
    """

    sample_xml = """    <mlrs:MiriLRS>
        <mlrs:AcqTargetID>Same Target as Observation</mlrs:AcqTargetID>
        <mlrs:AcqFilter>F1000W</mlrs:AcqFilter>
        <mlrs:AcqReadoutPattern>FAST</mlrs:AcqReadoutPattern>
        <mlrs:AcqGroups>3</mlrs:AcqGroups>
        <mlrs:AcqEtcId></mlrs:AcqEtcId>
        <mlrs:Subarray>FULL</mlrs:Subarray>
        <mlrs:Exposures>1</mlrs:Exposures>
        <mlrs:ReadoutPattern>FAST</mlrs:ReadoutPattern>
        <mlrs:Groups>10</mlrs:Groups>
        <mlrs:Integrations>3</mlrs:Integrations>
        <mlrs:EtcId></mlrs:EtcId>
        <mlrs:DitherType>MAPPING</mlrs:DitherType>
        <mlrs:NumberOfSpectralSteps>1</mlrs:NumberOfSpectralSteps>
        <mlrs:SpectralStepOffset>0</mlrs:SpectralStepOffset>
        <mlrs:NumberOfSpatialSteps>45</mlrs:NumberOfSpatialSteps>
        <mlrs:SpatialStepOffset>0.11</mlrs:SpatialStepOffset>
    </mlrs:MiriLRS>"""

    tree = read_fake_xml(sample_xml)

    metadata = {'obs_id': '1', 'filename': 'apt/1042.aptx', 'instrument': 'MIRI'}

    obj = templates.MiriLRS(tree, metadata=metadata)
    output_dict = obj.getobs()

    expected_dict = [{'ima_frames': 10, 'ima_integrations': 3, 'obs_id': '1', 'filename': 'apt/1042.aptx', 'instrument': 'MIRI',
                     'exposures': 1, 'readDetect': 'FAST', 'template': 'MiriLRS', 'NDither': 45, 'subarray': 'FULL'}]

    assertListDictEqual(output_dict, expected_dict)


def test_MiriMRS():
    """
    Test of MiriMRS template

    Sample from APT STScI ID 1031 observation 1
    Sample from APT STScI ID 1023 observation 23

    """

    sample_xml = """    <mmrs:MiriMRS>
        <mmrs:AcqTargetID>Same Target as Observation</mmrs:AcqTargetID>
        <mmrs:AcqFilter>F560W</mmrs:AcqFilter>
        <mmrs:AcqReadoutPattern>FAST</mmrs:AcqReadoutPattern>
        <mmrs:AcqGroups>5</mmrs:AcqGroups>
        <mmrs:AcqEtcId>5777</mmrs:AcqEtcId>
        <mmrs:Detector>MRS</mmrs:Detector>
        <mmrs:Dithers>
            <mmrs:MrsDitherSpecification>
                <mmrs:DitherType>4-Point</mmrs:DitherType>
                <mmrs:OptimizedFor>ALL</mmrs:OptimizedFor>
                <mmrs:Direction>Negative</mmrs:Direction>
            </mmrs:MrsDitherSpecification>
        </mmrs:Dithers>
        <mmrs:SimultaneousImaging>NO</mmrs:SimultaneousImaging>
        <mmrs:Subarray>FULL</mmrs:Subarray>
        <mmrs:PrimaryChannel>ALL</mmrs:PrimaryChannel>
        <mmrs:ExposureList>
            <mmrs:Exposure>
                <mmrs:Exposures>1</mmrs:Exposures>
                <mmrs:Wavelength>SHORT(A)</mmrs:Wavelength>
                <mmrs:Dither>Dither 1</mmrs:Dither>
                <mmrs:ReadoutPatternLong>FAST</mmrs:ReadoutPatternLong>
                <mmrs:GroupsLong>6</mmrs:GroupsLong>
                <mmrs:GroupsShort>6</mmrs:GroupsShort>
                <mmrs:IntegrationsLong>1</mmrs:IntegrationsLong>
                <mmrs:IntegrationsShort>1</mmrs:IntegrationsShort>
                <mmrs:EtcIdLong></mmrs:EtcIdLong>
                <mmrs:EtcIdShort></mmrs:EtcIdShort>
            </mmrs:Exposure>
            <mmrs:Exposure>
                <mmrs:Exposures>1</mmrs:Exposures>
                <mmrs:Wavelength>MEDIUM(B)</mmrs:Wavelength>
                <mmrs:Dither>Dither 1</mmrs:Dither>
                <mmrs:ReadoutPatternLong>FAST</mmrs:ReadoutPatternLong>
                <mmrs:GroupsLong>6</mmrs:GroupsLong>
                <mmrs:GroupsShort>6</mmrs:GroupsShort>
                <mmrs:IntegrationsLong>1</mmrs:IntegrationsLong>
                <mmrs:IntegrationsShort>1</mmrs:IntegrationsShort>
                <mmrs:EtcIdLong></mmrs:EtcIdLong>
                <mmrs:EtcIdShort></mmrs:EtcIdShort>
            </mmrs:Exposure>
            <mmrs:Exposure>
                <mmrs:Exposures>1</mmrs:Exposures>
                <mmrs:Wavelength>LONG(C)</mmrs:Wavelength>
                <mmrs:Dither>Dither 1</mmrs:Dither>
                <mmrs:ReadoutPatternLong>FAST</mmrs:ReadoutPatternLong>
                <mmrs:GroupsLong>6</mmrs:GroupsLong>
                <mmrs:GroupsShort>6</mmrs:GroupsShort>
                <mmrs:IntegrationsLong>1</mmrs:IntegrationsLong>
                <mmrs:IntegrationsShort>1</mmrs:IntegrationsShort>
                <mmrs:EtcIdLong></mmrs:EtcIdLong>
                <mmrs:EtcIdShort></mmrs:EtcIdShort>
            </mmrs:Exposure>
        </mmrs:ExposureList>
    </mmrs:MiriMRS>"""

    tree = read_fake_xml(sample_xml)

    metadata = {'obs_id': '1', 'filename': 'apt/1031.aptx', 'instrument': 'MIRI'}

    obj = templates.MiriMRS(tree, metadata=metadata)
    output_dict = obj.getobs()

    expected_dict = [{'LW_integrations': 1, 'disperser': 'SHORT', 'SW_frames': 6, 'obs_id': '1', 'filename': 'apt/1031.aptx',
                     'SW_integrations': 1, 'subarray': 'BOTH', 'simultaneous_imaging': 'NO', 'instrument': 'MIRI',
                     'detector': 'MRS', 'exposures': 1, 'LW_frames': 6, 'readDetect': 'FAST',
                     'template': 'MiriMRS', 'NDither': 4},
                     {'LW_integrations': 1, 'disperser': 'MEDIUM', 'SW_frames': 6, 'obs_id': '1',
                      'filename': 'apt/1031.aptx',
                      'SW_integrations': 1, 'subarray': 'BOTH', 'simultaneous_imaging': 'NO', 'instrument': 'MIRI',
                      'detector': 'MRS', 'exposures': 1, 'LW_frames': 6, 'readDetect': 'FAST',
                      'template': 'MiriMRS', 'NDither': 4},
                     {'LW_integrations': 1, 'disperser': 'LONG', 'SW_frames': 6, 'obs_id': '1',
                      'filename': 'apt/1031.aptx',
                      'SW_integrations': 1, 'subarray': 'BOTH', 'simultaneous_imaging': 'NO', 'instrument': 'MIRI',
                      'detector': 'MRS', 'exposures': 1, 'LW_frames': 6, 'readDetect': 'FAST',
                      'template': 'MiriMRS', 'NDither': 4}
                     ]

    assertListDictEqual(output_dict, expected_dict)


def test_MiriMRS_parallel():
    """
    Test of MiriMRS template with a parallel observation

    Sample from APT STScI ID 1024 observation 1
    """

    sample_xml = """    <mmrs:MiriMRS>
        <mmrs:AcqTargetID>Same Target as Observation</mmrs:AcqTargetID>
        <mmrs:AcqFilter>F1500W</mmrs:AcqFilter>
        <mmrs:AcqReadoutPattern>FAST</mmrs:AcqReadoutPattern>
        <mmrs:AcqGroups>5</mmrs:AcqGroups>
        <mmrs:AcqEtcId>17925</mmrs:AcqEtcId>
        <mmrs:Detector>ALL</mmrs:Detector>
        <mmrs:Dithers>
            <mmrs:MrsDitherSpecification>
                <mmrs:DitherType>4-Point</mmrs:DitherType>
                <mmrs:OptimizedFor>Channel 1</mmrs:OptimizedFor>
                <mmrs:Direction>Positive</mmrs:Direction>
            </mmrs:MrsDitherSpecification>
            <mmrs:MrsDitherSpecification>
                <mmrs:DitherType>4-Point</mmrs:DitherType>
                <mmrs:OptimizedFor>Channel 1</mmrs:OptimizedFor>
                <mmrs:Direction>Negative</mmrs:Direction>
            </mmrs:MrsDitherSpecification>
        </mmrs:Dithers>
        <mmrs:SimultaneousImaging>YES</mmrs:SimultaneousImaging>
        <mmrs:Subarray>FULL</mmrs:Subarray>
        <mmrs:PrimaryChannel>CHANNEL1</mmrs:PrimaryChannel>
        <mmrs:ExposureList>
            <mmrs:Exposure>
                <mmrs:Exposures>1</mmrs:Exposures>
                <mmrs:Filter>F770W</mmrs:Filter>
                <mmrs:Wavelength>SHORT(A)</mmrs:Wavelength>
                <mmrs:Dither>Dither 1</mmrs:Dither>
                <mmrs:ReadoutPattern>FAST</mmrs:ReadoutPattern>
                <mmrs:ReadoutPatternLong>FAST</mmrs:ReadoutPatternLong>
                <mmrs:ReadoutPatternShort>FAST</mmrs:ReadoutPatternShort>
                <mmrs:Groups>10</mmrs:Groups>
                <mmrs:GroupsLong>10</mmrs:GroupsLong>
                <mmrs:GroupsShort>10</mmrs:GroupsShort>
                <mmrs:Integrations>10</mmrs:Integrations>
                <mmrs:IntegrationsLong>10</mmrs:IntegrationsLong>
                <mmrs:IntegrationsShort>10</mmrs:IntegrationsShort>
                <mmrs:EtcId>17925</mmrs:EtcId>
                <mmrs:EtcIdLong></mmrs:EtcIdLong>
                <mmrs:EtcIdShort></mmrs:EtcIdShort>
            </mmrs:Exposure>
            <mmrs:Exposure>
                <mmrs:Exposures>1</mmrs:Exposures>
                <mmrs:Filter>F770W</mmrs:Filter>
                <mmrs:Wavelength>SHORT(A)</mmrs:Wavelength>
                <mmrs:Dither>Dither 2</mmrs:Dither>
                <mmrs:ReadoutPattern>FAST</mmrs:ReadoutPattern>
                <mmrs:ReadoutPatternLong>FAST</mmrs:ReadoutPatternLong>
                <mmrs:ReadoutPatternShort>FAST</mmrs:ReadoutPatternShort>
                <mmrs:Groups>10</mmrs:Groups>
                <mmrs:GroupsLong>10</mmrs:GroupsLong>
                <mmrs:GroupsShort>10</mmrs:GroupsShort>
                <mmrs:Integrations>10</mmrs:Integrations>
                <mmrs:IntegrationsLong>10</mmrs:IntegrationsLong>
                <mmrs:IntegrationsShort>10</mmrs:IntegrationsShort>
                <mmrs:EtcId></mmrs:EtcId>
                <mmrs:EtcIdLong></mmrs:EtcIdLong>
                <mmrs:EtcIdShort></mmrs:EtcIdShort>
            </mmrs:Exposure>
        </mmrs:ExposureList>
    </mmrs:MiriMRS>"""

    tree = read_fake_xml(sample_xml)

    metadata = {'obs_id': '1', 'filename': 'apt/1024.aptx', 'instrument': 'MIRI'}

    obj = templates.MiriMRS(tree, metadata=metadata)
    output_dict = obj.getobs()

    expected_dict = [
                    {'LW_frames': 10, 'LW_integrations': 10, 'NDither': 4, 'SW_frames': 10,
                     'SW_integrations': 10, 'detector': 'ALL', 'disperser': 'SHORT',
                    'exposures': 1, 'filename': 'apt/1024.aptx', 'filter': 'F770W',
                    'ima_frames': 10, 'ima_integrations': 10, 'instrument': 'MIRI',
                    'obs_id': '1', 'primary_channel': 'CHANNEL1', 'readDetect': 'FAST',
                    'simultaneous_imaging': 'YES', 'subarray': 'FULL', 'template': 'MiriMRS'},
                     {'LW_frames': 10, 'LW_integrations': 10, 'NDither': 4, 'SW_frames': 10,
                     'SW_integrations': 10, 'detector': 'ALL', 'disperser': 'SHORT',
                     'exposures': 1, 'filename': 'apt/1024.aptx', 'filter': 'F770W',
                     'ima_frames': 10, 'ima_integrations': 10, 'instrument': 'MIRI',
                     'obs_id': '1', 'primary_channel': 'CHANNEL1', 'readDetect': 'FAST',
                     'simultaneous_imaging': 'YES', 'subarray': 'FULL', 'template': 'MiriMRS'}
                     ]

    assertListDictEqual(output_dict, expected_dict)


def test_MiriMRSCrossGratingEngineering():
    """
    Test of MiriMRSCrossGratingEngineering template

    Sample from APT STScI ID 1050 observation 2
    """

    sample_xml = """    <mmrscge:MiriMRSCrossGratingEngineering>
        <mmrscge:ExposureList>
            <mmrscge:Exposure>
                <mmrscge:IsLampOn>false</mmrscge:IsLampOn>
                <mmrscge:Exposures>1</mmrscge:Exposures>
                <mmrscge:IntegrationsLong>10</mmrscge:IntegrationsLong>
                <mmrscge:IntegrationsShort>20</mmrscge:IntegrationsShort>
                <mmrscge:GroupsShort>5</mmrscge:GroupsShort>
                <mmrscge:GroupsLong>10</mmrscge:GroupsLong>
                <mmrscge:EtcId></mmrscge:EtcId>
            </mmrscge:Exposure>
            <mmrscge:Exposure>
                <mmrscge:IsLampOn>false</mmrscge:IsLampOn>
                <mmrscge:Exposures>1</mmrscge:Exposures>
                <mmrscge:IntegrationsLong>1</mmrscge:IntegrationsLong>
                <mmrscge:IntegrationsShort>1</mmrscge:IntegrationsShort>
                <mmrscge:GroupsShort>2</mmrscge:GroupsShort>
                <mmrscge:GroupsLong>2</mmrscge:GroupsLong>
                <mmrscge:EtcId></mmrscge:EtcId>
            </mmrscge:Exposure>
            <mmrscge:Exposure>
                <mmrscge:IsLampOn>false</mmrscge:IsLampOn>
                <mmrscge:Exposures>1</mmrscge:Exposures>
                <mmrscge:IntegrationsLong>1</mmrscge:IntegrationsLong>
                <mmrscge:IntegrationsShort>1</mmrscge:IntegrationsShort>
                <mmrscge:GroupsShort>2</mmrscge:GroupsShort>
                <mmrscge:GroupsLong>2</mmrscge:GroupsLong>
                <mmrscge:EtcId></mmrscge:EtcId>
            </mmrscge:Exposure>
            <mmrscge:Exposure>
                <mmrscge:IsLampOn>false</mmrscge:IsLampOn>
                <mmrscge:Exposures>1</mmrscge:Exposures>
                <mmrscge:IntegrationsLong>10</mmrscge:IntegrationsLong>
                <mmrscge:IntegrationsShort>20</mmrscge:IntegrationsShort>
                <mmrscge:GroupsShort>5</mmrscge:GroupsShort>
                <mmrscge:GroupsLong>10</mmrscge:GroupsLong>
                <mmrscge:EtcId></mmrscge:EtcId>
            </mmrscge:Exposure>
            <mmrscge:Exposure>
                <mmrscge:IsLampOn>false</mmrscge:IsLampOn>
                <mmrscge:Exposures>1</mmrscge:Exposures>
                <mmrscge:IntegrationsLong>1</mmrscge:IntegrationsLong>
                <mmrscge:IntegrationsShort>1</mmrscge:IntegrationsShort>
                <mmrscge:GroupsShort>2</mmrscge:GroupsShort>
                <mmrscge:GroupsLong>2</mmrscge:GroupsLong>
                <mmrscge:EtcId></mmrscge:EtcId>
            </mmrscge:Exposure>
            <mmrscge:Exposure>
                <mmrscge:IsLampOn>false</mmrscge:IsLampOn>
                <mmrscge:Exposures>1</mmrscge:Exposures>
                <mmrscge:IntegrationsLong>1</mmrscge:IntegrationsLong>
                <mmrscge:IntegrationsShort>1</mmrscge:IntegrationsShort>
                <mmrscge:GroupsShort>2</mmrscge:GroupsShort>
                <mmrscge:GroupsLong>2</mmrscge:GroupsLong>
                <mmrscge:EtcId></mmrscge:EtcId>
            </mmrscge:Exposure>
            <mmrscge:Exposure>
                <mmrscge:IsLampOn>false</mmrscge:IsLampOn>
                <mmrscge:Exposures>1</mmrscge:Exposures>
                <mmrscge:IntegrationsLong>10</mmrscge:IntegrationsLong>
                <mmrscge:IntegrationsShort>20</mmrscge:IntegrationsShort>
                <mmrscge:GroupsShort>5</mmrscge:GroupsShort>
                <mmrscge:GroupsLong>10</mmrscge:GroupsLong>
                <mmrscge:EtcId></mmrscge:EtcId>
            </mmrscge:Exposure>
            <mmrscge:Exposure>
                <mmrscge:IsLampOn>false</mmrscge:IsLampOn>
                <mmrscge:Exposures>1</mmrscge:Exposures>
                <mmrscge:IntegrationsLong>1</mmrscge:IntegrationsLong>
                <mmrscge:IntegrationsShort>1</mmrscge:IntegrationsShort>
                <mmrscge:GroupsShort>2</mmrscge:GroupsShort>
                <mmrscge:GroupsLong>2</mmrscge:GroupsLong>
                <mmrscge:EtcId></mmrscge:EtcId>
            </mmrscge:Exposure>
            <mmrscge:Exposure>
                <mmrscge:IsLampOn>false</mmrscge:IsLampOn>
                <mmrscge:Exposures>1</mmrscge:Exposures>
                <mmrscge:IntegrationsLong>1</mmrscge:IntegrationsLong>
                <mmrscge:IntegrationsShort>1</mmrscge:IntegrationsShort>
                <mmrscge:GroupsShort>2</mmrscge:GroupsShort>
                <mmrscge:GroupsLong>2</mmrscge:GroupsLong>
                <mmrscge:EtcId></mmrscge:EtcId>
            </mmrscge:Exposure>
        </mmrscge:ExposureList>
        <mmrscge:LampUse>OFF ONLY</mmrscge:LampUse>
    </mmrscge:MiriMRSCrossGratingEngineering>"""

    tree = read_fake_xml(sample_xml)

    metadata = {'obs_id': '2', 'filename': 'apt/1050.aptx', 'instrument': 'MIRI'}

    obj = templates.MiriMRSCrossGratingEngineering(tree, metadata=metadata)
    output_dict = obj.getobs()

    expected_dict = [{'LW_integrations': 10, 'LW_frames': 10, 'obs_id': '2', 'filename': 'apt/1050.aptx',
        'SW_integrations': 20, 'SW_frames': 5, 'subarray': 'FULL', 'detector': 'BOTH', 'exposures': 1,
        'template': 'MiriMRSCrossGratingEngineering', 'instrument': 'MIRI', 'NDither': 1},
        {'LW_integrations': 1, 'LW_frames': 2, 'obs_id': '2', 'filename': 'apt/1050.aptx',
        'SW_integrations': 1, 'SW_frames': 2, 'subarray': 'FULL', 'detector': 'BOTH', 'exposures': 1,
        'template': 'MiriMRSCrossGratingEngineering', 'instrument': 'MIRI', 'NDither': 1},
        {'LW_integrations': 1, 'LW_frames': 2, 'obs_id': '2', 'filename': 'apt/1050.aptx',
        'SW_integrations': 1, 'SW_frames': 2, 'subarray': 'FULL', 'detector': 'BOTH', 'exposures': 1,
        'template': 'MiriMRSCrossGratingEngineering', 'instrument': 'MIRI', 'NDither': 1},
        {'LW_integrations': 10, 'LW_frames': 10, 'obs_id': '2', 'filename': 'apt/1050.aptx',
        'SW_integrations': 20, 'SW_frames': 5, 'subarray': 'FULL', 'detector': 'BOTH', 'exposures': 1,
        'template': 'MiriMRSCrossGratingEngineering', 'instrument': 'MIRI', 'NDither': 1},
        {'LW_integrations': 1, 'LW_frames': 2, 'obs_id': '2', 'filename': 'apt/1050.aptx',
        'SW_integrations': 1, 'SW_frames': 2, 'subarray': 'FULL', 'detector': 'BOTH', 'exposures': 1,
        'template': 'MiriMRSCrossGratingEngineering', 'instrument': 'MIRI', 'NDither': 1},
        {'LW_integrations': 1, 'LW_frames': 2, 'obs_id': '2', 'filename': 'apt/1050.aptx',
        'SW_integrations': 1, 'SW_frames': 2, 'subarray': 'FULL', 'detector': 'BOTH', 'exposures': 1,
        'template': 'MiriMRSCrossGratingEngineering', 'instrument': 'MIRI', 'NDither': 1},
        {'LW_integrations': 10, 'LW_frames': 10, 'obs_id': '2', 'filename': 'apt/1050.aptx',
        'SW_integrations': 20, 'SW_frames': 5, 'subarray': 'FULL', 'detector': 'BOTH', 'exposures': 1,
        'template': 'MiriMRSCrossGratingEngineering', 'instrument': 'MIRI', 'NDither': 1},
        {'LW_integrations': 1, 'LW_frames': 2, 'obs_id': '2', 'filename': 'apt/1050.aptx',
        'SW_integrations': 1, 'SW_frames': 2, 'subarray': 'FULL', 'detector': 'BOTH', 'exposures': 1,
        'template': 'MiriMRSCrossGratingEngineering', 'instrument': 'MIRI', 'NDither': 1},
        {'LW_integrations': 1, 'LW_frames': 2, 'obs_id': '2', 'filename': 'apt/1050.aptx',
        'SW_integrations': 1, 'SW_frames': 2, 'subarray': 'FULL', 'detector': 'BOTH', 'exposures': 1,
        'template': 'MiriMRSCrossGratingEngineering', 'instrument': 'MIRI', 'NDither': 1}]

    assertListDictEqual(output_dict, expected_dict)


if __name__ == '__main__':
    import sys
    import inspect

    # List all functions starting with test_ to launch them automatically
    testfunctions = [obj for name, obj in inspect.getmembers(sys.modules[__name__])
                     if (inspect.isfunction(obj) and
                         name.startswith('test_'))]

    for func in testfunctions:
        print("Running {}".format(func.__name__))
        func()