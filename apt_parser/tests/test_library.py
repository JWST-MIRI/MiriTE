import pytest
import numpy.testing as npt
from miri.apt_parser import library
from miri.apt_parser.utils import assertDictEqual, read_test_xml

testtarget = [("""        <apt:Target xsi:type="FixedTargetType" AutoGenerated="false" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">
            <apt:Number>1</apt:Number>
            <apt:TargetName>LMC-STAR-A</apt:TargetName>
            <apt:TargetArchiveName>2MASS J05215224-6930510</apt:TargetArchiveName>
            <apt:TargetID>LMC-STAR-A</apt:TargetID>
            <apt:Comments>Coordinates from Anderson's JWST astrometric catalogue</apt:Comments>
            <apt:RAProperMotion></apt:RAProperMotion>
            <apt:DecProperMotion></apt:DecProperMotion>
            <apt:RAProperMotionUnits>mas/yr</apt:RAProperMotionUnits>
            <apt:DecProperMotionUnits>mas/yr</apt:DecProperMotionUnits>
            <apt:Epoch>2000</apt:Epoch>
            <apt:Extended>NO</apt:Extended>
            <apt:Category>Calibration</apt:Category>
            <apt:Keywords>Aperture location</apt:Keywords>
            <apt:Keywords>Astrometric</apt:Keywords>
            <apt:Keywords>External flat field</apt:Keywords>
            <apt:Keywords>Telescope/sky background</apt:Keywords>
            <apt:EquatorialCoordinates Value="05 21 52.2463 -69 30 51.05">
                <apt:RAUncertainty>
                    <apt:Value></apt:Value>
                    <apt:Units>Arcsec</apt:Units>
                </apt:RAUncertainty>
                <apt:DecUncertainty>
                    <apt:Value></apt:Value>
                    <apt:Units>Arcsec</apt:Units>
                </apt:DecUncertainty>
            </apt:EquatorialCoordinates>
            <apt:BackgroundTargetReq>false</apt:BackgroundTargetReq>
            <apt:TargetConfirmationRun>true</apt:TargetConfirmationRun>
        </apt:Target>""", # APT ID 1024, target #1
{'name': 'LMC-STAR-A', 'number': 1, 'background': False}),
              ("""        <apt:Target xsi:type="FixedTargetType" AutoGenerated="false" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">
            <apt:Number>2</apt:Number>
            <apt:TargetName>LMC-STAR-BACKUP</apt:TargetName>
            <apt:TargetArchiveName>2MASS J05220207-6930388</apt:TargetArchiveName>
            <apt:TargetID>LMC-STAR-BACKUP</apt:TargetID>
            <apt:Comments>Coordinates from Anderson's JWST astrometric catalogue</apt:Comments>
            <apt:RAProperMotion>25.0</apt:RAProperMotion>
            <apt:DecProperMotion>-43.5</apt:DecProperMotion>
            <apt:RAProperMotionUnits>mas/yr</apt:RAProperMotionUnits>
            <apt:DecProperMotionUnits>mas/yr</apt:DecProperMotionUnits>
            <apt:Epoch>2000</apt:Epoch>
            <apt:Extended>NO</apt:Extended>
            <apt:Category>Calibration</apt:Category>
            <apt:Keywords>Aperture location</apt:Keywords>
            <apt:Keywords>Astrometric</apt:Keywords>
            <apt:Keywords>External flat field</apt:Keywords>
            <apt:Keywords>Telescope/sky background</apt:Keywords>
            <apt:EquatorialCoordinates Value="05 22 2.0669 -69 30 38.74">
                <apt:RAUncertainty>
                    <apt:Value></apt:Value>
                    <apt:Units>Arcsec</apt:Units>
                </apt:RAUncertainty>
                <apt:DecUncertainty>
                    <apt:Value></apt:Value>
                    <apt:Units>Arcsec</apt:Units>
                </apt:DecUncertainty>
            </apt:EquatorialCoordinates>
            <apt:BackgroundTargetReq>false</apt:BackgroundTargetReq>
            <apt:TargetConfirmationRun>true</apt:TargetConfirmationRun>
        </apt:Target>""", # APT 1024, target #2
               {'name': 'LMC-STAR-BACKUP', 'number': 2, 'background': False}),
              ("""<apt:Target xsi:type="FixedTargetType" AutoGenerated="false" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">
            <apt:Number>1</apt:Number>
            <apt:TargetName>GJ-1214</apt:TargetName>
            <apt:TargetArchiveName>GJ 1214</apt:TargetArchiveName>
            <apt:TargetID>GJ-1214</apt:TargetID>
            <apt:Comments>This object was generated by the targetselector and retrieved from the SIMBAD database.</apt:Comments>
            <apt:RAProperMotion></apt:RAProperMotion>
            <apt:DecProperMotion></apt:DecProperMotion>
            <apt:RAProperMotionUnits></apt:RAProperMotionUnits>
            <apt:DecProperMotionUnits></apt:DecProperMotionUnits>
            <apt:Extended>Unknown</apt:Extended>
            <apt:Category>Star</apt:Category>
            <apt:EquatorialCoordinates Value="17 15 18.9400 +04 57 49.70">
                <apt:RAUncertainty>
                    <apt:Value></apt:Value>
                    <apt:Units>Arcsec</apt:Units>
                </apt:RAUncertainty>
                <apt:DecUncertainty>
                    <apt:Value></apt:Value>
                    <apt:Units>Arcsec</apt:Units>
                </apt:DecUncertainty>
            </apt:EquatorialCoordinates>
            <apt:BackgroundTargetReq>true</apt:BackgroundTargetReq>
            <apt:BackgroundTargets>2 GJ-1214-BK</apt:BackgroundTargets>
            <apt:BackgroundTargets>3 GJ-1214-BK-COPY</apt:BackgroundTargets>
            <apt:TargetConfirmationRun>false</apt:TargetConfirmationRun>
        </apt:Target>""",
               {'name': 'GJ-1214', 'number': 1, 'background': True, "bkg_targets":["2 GJ-1214-BK", "3 GJ-1214-BK-COPY"]}),
              ]

def test_get_simplified_tag():
    input = "{http://www.stsci.edu/JWST/APT/Template/MiriLRS}MiriLRS"
    expected = "MiriLRS"

    output = library.get_simplified_tag(input)

    assert output == expected


@pytest.mark.parametrize("sample_xml,expected", testtarget)
def test_get_target(sample_xml, expected):
    """
    From Proposal ID 1024 (1st target)

    Note that for each tag, you must add "apt:" in front of it (both start and end)
    for the test to mimick what happens with a real XML file. Don't ask me why.

    :return:
    """

    tree = read_test_xml(sample_xml)

    target = library.get_target(tree)

    assertDictEqual(target, expected)


(ref_ram, ref_time, dummy) = library.get_prediction({'integrations':1, 'exposures':1, 'frames':5, "subarray": 'FULL', "NDither":1})

test_prediction = [
    ({'integrations':1, 'exposures':1, 'frames':5, "subarray": 'FULL', "NDither":4}, (ref_ram, ref_time, 4)),
    ({'integrations':1, 'exposures':1, 'frames':5, "subarray": 'BOTH', "NDither":4}, (ref_ram, ref_time, 8)),
    ({'integrations':1, 'exposures':3, 'frames':5, "subarray": 'FULL', "NDither":2}, (ref_ram, ref_time, 6)),
    ({'integrations':1, 'exposures':1, 'frames':5, "subarray": 'BRIGHTSKY', "NDither":1}, (ref_ram*0.24806201550387597, ref_time*0.24806201550387597, 1)),
]

@pytest.mark.parametrize("input,expected", test_prediction)
def test_get_prediction(input, expected):
    (ref_ram, ref_time, ref_exps) = expected

    (ram, time, exps) = library.get_prediction(input)

    npt.assert_almost_equal(ram, ref_ram, err_msg="RAM prediction different than expected")

    npt.assert_almost_equal(time, ref_time, err_msg="Time prediction different than expected")

    assert exps == ref_exps, "In get_prediction, Number of exposures ({}) different than expected ({}).".format(exps, ref_exps)