import logging
import re
import numpy as np

import xml.etree.ElementTree as ET
from .constants import ns

LOG = logging.getLogger(__name__)


def assert_list_dict_equal(input_list, ref_list):
    """
    Assert if All the dicts in both lists are equal (comparing indexes to indexes)

    :param list input_list: input list of dictionnaries
    :param list ref_list: reference list of dictionnaries

    """

    assert len(input_list) == len(ref_list), "List have different sizes"

    for (d1, d2) in zip(input_list, ref_list):
        assert_dict_equal(d1, d2)


def assert_dict_equal(input_dict, dict_ref):
    """
    Assert if Dicts are equal. If not, display the differences

    :param dict input_dict: Input dictionary
    :param dict dict_ref: Reference dictionary

    """

    (is_equal, msg) = compare_dict(input_dict, dict_ref)

    assert is_equal, msg


def compare_dict(input_dict, dict_ref, msg=None, prefix=None):
    r"""
    Compare 2 dicts. Will return an error message explaining the differences

    :param dict input_dict: Input dictionary
    :param dict dict_ref: Reference dictionary

    /!\ Warning: All other parameters are internal (for recursivity) and must NOT be used

    :return: if dict are equal and error message associated with it
    :rtype: (bool, msg)
    """

    is_equal = True

    if not msg:
        msg = ""

    keys1 = set(input_dict.keys())
    keys2 = set(dict_ref.keys())

    d1_prefix = "input_dict"
    d2_prefix = "dict_ref"
    if prefix:
        d1_prefix += prefix
        d2_prefix += prefix

    common_keys = keys1.intersection(keys2)

    # Keys present in keys1 not present in keys2
    new_keys1 = keys1.difference(keys2)
    new_keys1 = list(new_keys1)
    new_keys1.sort()
    if len(new_keys1) != 0:
        is_equal = False
        msg += "Keys exclusive to {}:\n".format(d1_prefix)
        for key in new_keys1:
            msg += "\t{}[{}] = {}\n".format(d1_prefix, key, input_dict[key])

    # Keys present in keys2 not present in keys1
    new_keys2 = keys2.difference(keys1)
    new_keys2 = list(new_keys2)
    new_keys2.sort()
    if len(new_keys2) != 0:
        is_equal = False
        msg += "Keys exclusive to {}:\n".format(d2_prefix)
        for key in new_keys2:
            msg += "\t{}[{}] = {}\n".format(d2_prefix, key, dict_ref[key])

    common_keys = list(common_keys)
    common_keys.sort()
    for key in common_keys:
        value1 = input_dict[key]
        value2 = dict_ref[key]
        if isinstance(value1, dict):
            new_prefix = prefix if prefix else ""
            new_prefix += "[{}]".format(key)
            (value_equal, tmp_msg) = compare_dict(value1, value2, prefix=new_prefix)
            if not value_equal:
                is_equal = False
            msg += tmp_msg
        else:
            if value1 != value2:
                is_equal = False
                msg += "Difference for:\n"
                msg += "\t{}[{}] = {}\n".format(d1_prefix, key, value1)
                msg += "\t{}[{}] = {}\n".format(d2_prefix, key, value2)

    return is_equal, msg


def read_fake_xml(txt):
    """
    Read an XML string, adding the generic namespace for the JWST APT in it and any other parts needed

    This function is designed solely for unit. testing since we need to parse only a fraction of an
    XML and not the entire one for unit testing

    The input string just have to be an Element (and ONLY one) and its sub-element if needed.

    We need that because the namespace needs to be added; and we want a simple way to copy paste code from an APT file
    without anything to add on our side (minimum maintenance effort)

    :param str txt: A string containing a random XML sub-part you want to parse and read.
                    This string MUST NOT have any white line before or after

    :return: XML ElementTree for that string
    :rtype:
    """
    namespace_list = []
    for (k, v) in ns.items():
        namespace_list.append('xmlns:{}="{}"'.format(k, v))

    xml_txt = '<?xml version="1.0" encoding="UTF-8"?>\n'
    xml_txt += "<JwstProposal {}>\n".format(" ".join(namespace_list))
    xml_txt += txt
    xml_txt += "\n</JwstProposal>"

    tmp = ET.fromstring(xml_txt)
    tree = ET.ElementTree(tmp)
    root = tree.getroot()

    # Get the actual XML code, on the condition the user provided only one item in TXT
    el = root.find("./")

    return el


def get_namespace(xml_text):
    """
    Given an XML text, and assuming there is only one namespace in it for all the tags, return that short namespace
    as a string


    :param str xml_text: xml sample text
    :return: namespace as a string (for instance mmrs)
    :rtype: str
    """

    pattern = re.compile(r"<([^\/:]+):")

    match = pattern.findall(xml_text)

    if match:
        return match[0]
    else:
        LOG.error("No namespace could be found in {}".format(xml_text))
        raise ValueError("No namespace could be found in {}".format(xml_text))


def strtime(sec):
    """
    Format a time interval in seconds.
    Return a string like:
    42 days 15h25m54s

    Note: The number of days is omitted if it's 0

    Return "N/A" if not a float/int or NaN

    :param float sec: Time interval in seconds
    :return str: Formatted string to display this time nicely
    """

    if not isinstance(sec, (int, float)) or np.isnan(sec):
        return "N/A"

    (d, r) = divmod(sec, 86400)
    (h, r) = divmod(r, 3600)
    (m, s) = divmod(r, 60)

    if d:
        string = "{:.0f} days ".format(d)
    else:
        string = ""

    string += "{:02.0f}h{:02.0f}m{:02.0f}s".format(h, m, s)

    return string
