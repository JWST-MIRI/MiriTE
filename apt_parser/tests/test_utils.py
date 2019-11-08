import pytest
from miri.apt_parser import utils

def test_get_namespace():
    sample_xml = """    <mmrs:MrsDitherSpecification>
        <mmrs:DitherType>4-Point</mmrs:DitherType>
        <mmrs:OptimizedFor>Channel 1</mmrs:OptimizedFor>
        <mmrs:Direction>Positive</mmrs:Direction>
    </mmrs:MrsDitherSpecification>"""

    namespace = utils.get_namespace(sample_xml)

    expected_namespace = "mmrs"

    assert namespace == expected_namespace, "Namespace: Expected {}, got {}".format(expected_namespace, namespace)

    sample_xml = "   </mmrs:MrsDitherSpecification>"
    with pytest.raises(ValueError):
        namespace = utils.get_namespace(sample_xml)
