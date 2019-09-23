import xml.etree.ElementTree as ET
import xml
import logging
import numpy as np
import os
import zipfile
import collections
import re

LOG = logging.getLogger('parse_apt.library')

from . import templates as tp
from . import constants as c


def get_simplified_tag(tag):
    """

    For an input tag like:
    {http://www.stsci.edu/JWST/APT/Template/MiriLRS}MiriLRS

    return:
    "MiriLRS"

    :param tag:
    :type tag:
    :return:
    :rtype:
    """
    regexp = r"^\{.*\}(.+)$"

    pattern_obj = re.compile(regexp)

    match = pattern_obj.match(tag)

    if match:
        return match.group(1)
    else:
        raise ValueError("Regexp {} don't match in '{}'".format(regexp, tag))

def get_element_tree(filename):
    """
    return xml.etree.ElementTree from the input file (.xml or .aptx)

    :param str filename: XML filename (either .xml or .aptx)
    :return: XML tree object
    :rtype: xml.etree.ElementTree
    """

    basename = os.path.basename(filename)
    (bare_name, ext) = os.path.splitext(basename)

    if ext == ".xml":
        try:
            tree = ET.parse(filename)
            root = tree.getroot()
        except xml.etree.ElementTree.ParseError:
            LOG.error("XML parsing failed for {}, skipping.".format(filename))
            return None
    elif ext == ".aptx":
        zf = zipfile.ZipFile(filename, 'r')
        f = zf.read("{}.xml".format(bare_name))
        root = ET.fromstring(f)
    else:
        LOG.error("Unknown file extension: {} ({})".format(ext, filename))

    return root


def get_target(template):
    """

    :param template: MiriImaging observation XML object
    :type template: <class 'xml.etree.ElementTree.Element'>

    :return:
    """

    target_number = int(template.find("apt:Number", c.ns).text)
    target_name = template.find("apt:TargetName", c.ns).text

    bool_txt = template.find("apt:BackgroundTargetReq", c.ns).text

    if bool_txt == "true":
        background = True
    elif bool_txt == "false":
        background = False
    else:
        raise ValueError("Unknown value for bool_txt: {}".format(bool_txt))

    output_dict = {'name': target_name, 'number': target_number, 'background': background}

    if background:
        bkg_targets = [item.text for item in template.findall("apt:BackgroundTargets", c.ns)]
        output_dict["bkg_targets"] = bkg_targets

    return output_dict

def parse_apt(filename):
    """
    Parse an APT file given its filename (can be either .aptx or .apt/XML)

    The output will be a list of individual simulations for that APT file.


    :param str filename:
    :return: properties of the APT
    :rtype: dict
    """

    # List of tuples (exps, ints, groups)
    simulation_list = []

    root = get_element_tree(filename)

    skip_obs = ["MiriAnneal"]


    # Select all targets
    xml_targets = root.findall('.//apt:Targets/apt:Target', c.ns)

    targets = {}
    targets["NONE"] = {'name': None, "background": None}  # Default target for template where it's irrelevant
    for xml_target in xml_targets:
        target = get_target(xml_target)
        targets["{} {}".format(target["number"], target["name"])] = target

    # Select all observation regardless of their position from root
    observations = root.findall('.//apt:ObservationGroup/apt:Observation', c.ns)

    #miri_imaging_list = root.findall('.//apt:Observation//mi:MiriImaging', parser.ns)



    for obs in observations:
        obs_id = obs.find("apt:Number", c.ns).text
        target_id = obs.find("apt:TargetID", c.ns).text
        target = targets[target_id]
        instrument = obs.find("apt:Instrument", c.ns).text

        # we take the first child of Template as we don't expect something else
        template = obs.find("apt:Template//", c.ns)


        #MIRIImaging
        tag = get_simplified_tag(template.tag)

        if tag in skip_obs:
            continue

        metadata = {
            "obs_id":obs_id,
            "instrument":instrument,
            "filename":filename,
            "target_id":target_id,
            "background": target["background"]

        }

        try:
            obj = tp.templates[tag](template, metadata=metadata)
            tmp_sim_list = obj.getobs()
            simulation_list.extend(tmp_sim_list)

        except KeyError:
            LOG.info("{}: Observation {}".format(filename, template.tag))




    return simulation_list

def get_prediction(sim):
    """
    Predict time in hours and memory in GB needed for the input simulation

    :param dict sim: Dictionnary containing metadata for the simulations.
                     (integrations, frames, subarray, exposures and NDither are needed)

    :return: (Ram in GB, time in hours, Number of exposures)
    :rtype: (float, float, int)
    """

    # Can be in a key named: integrations, ima_integrations, LW_integrations, SW_integrations
    keys = list(filter(lambda x: 'integrations' in x, sim.keys()))
    if len(keys) == 1:
        key_prefix = keys[0].rstrip("integrations")
    else:
        raise ValueError("Unable to find nb_integrations in {} (too many or no corresponding value)".format(sim))

    integration_key = "{}integrations".format(key_prefix)
    frame_key = "{}frames".format(key_prefix)

    # Values extracted from my benchmark, done with Pipeline 7.3 (July 2019)
    ram = 0.05 * sim[integration_key] * sim[frame_key] + 1.7
    time = (8.75 * sim[integration_key] * sim[frame_key] + 50.)/ 3600.  # In hours

    # Correct for detector size (compared to FULL ARRAY)
    corr_factor = c.SUBARRAY_PIX[sim["subarray"]] / (1032. * 1024.)
    ram *= corr_factor
    time *= corr_factor

    nb_exps = sim["exposures"] * sim["NDither"]

    if sim["subarray"] == "BOTH":
        nb_exps *= 2

    return ram, time, nb_exps

def count_datavolume(sim_dict):
    """
    Extract from the given input the amount of time and the memory you need to
    process each simulation through the JWST pipeline

    :param dict sim_dict: Each key represent a set of simulations (a CAR activity for instance)
                          each value is a list of simulation. Each simulation being a dict with detailled info
    :return: Return (mem, time) where mem and time are dictionnaries with the same keys as the input dict.
    :rtype: Memory is in GB, Time is in hours
    """
    mem_volume = {}  # Total memory required in GB
    time_volume = {}  # Pipeline estimated run time in s

    for (car, sim_list) in sim_dict.items():
        memory = []
        times = []
        for sim in sim_list:
            if "detector" in sim.keys():
                if sim["detector"] in ["IMAGER", "ALL"]:
                    tmp = {
                        "integrations":sim["ima_integrations"],
                        "frames":sim["ima_frames"],
                        "exposures":sim["exposures"],
                        "subarray":sim["subarray"],
                        "NDither":sim["NDither"],
                    }
                    (ram, time, nb_exps) = get_prediction(tmp)
                    memory.extend([ram] * nb_exps)  # For each exposure we have one identical file to analyse
                    times.extend([time] * nb_exps)  # For each exposure we have one identical file to analyse

                if sim["detector"] in ["ALL", "MRS"]:
                    tmp = {
                        "integrations": sim["LW_integrations"],
                        "frames": sim["LW_frames"],
                        "exposures": sim["exposures"],
                        "subarray": "FULL",
                        "NDither": sim["NDither"],
                    }
                    (ram, time, nb_exps) = get_prediction(tmp)
                    memory.extend([ram] * nb_exps)  # For each exposure we have one identical file to analyse
                    times.extend([time] * nb_exps)  # For each exposure we have one identical file to analyse

                    tmp = {
                        "integrations": sim["SW_integrations"],
                        "frames": sim["SW_frames"],
                        "exposures": sim["exposures"],
                        "subarray": "FULL",
                        "NDither": sim["NDither"],
                    }
                    (ram, time, nb_exps) = get_prediction(tmp)
                    memory.extend([ram] * nb_exps)  # For each exposure we have one identical file to analyse
                    times.extend([time] * nb_exps)  # For each exposure we have one identical file to analyse


            else:
                (ram, time, nb_exps) = get_prediction(sim)
                memory.extend([ram] * nb_exps)  # For each exposure we have one identical file to analyse
                times.extend([time] * nb_exps)  # For each exposure we have one identical file to analyse

        mem_volume[car] = np.array(memory)
        time_volume[car] = np.array(times)

    return mem_volume, time_volume


def init_log(log="parse_apt.log", stdout_loglevel="INFO", file_loglevel="DEBUG", extra_config=None):
    """
    Init logging configuration file. Must be used before any other package you could think of and that might use logging as well.
    If it doesn't work the way you expect, try inverting the imports to see if it changes.

    :param str log: filename where to store logs. By default "pipeline.log"
    :param str stdout_loglevel: log level for standard output (ERROR, WARNING, INFO, DEBUG)
    :param str file_loglevel: log level for log file (ERROR, WARNING, INFO, DEBUG)
    :param dict extra_config: [optional] Set of extra properties to be added to the dict_config for logging
    :return:
    :rtype:
    """

    import logging.config

    log_config = {
        "version": 1,
        "formatters":
            {
                "form01":
                    {
                        "format": "%(asctime)s %(levelname)-8s %(message)s",
                        "datefmt": "%H:%M:%S"
                    },
                "form02":
                    {
                        "format": "%(asctime)s [%(processName)s/%(name)s] %(levelname)s - %(message)s",
                        "datefmt": "%H:%M:%S"
                    },
            },
        "handlers":
            {
                "console":
                    {
                        "class": "logging.StreamHandler",
                        "formatter": "form01",
                        "level": stdout_loglevel,
                        "stream": "ext://sys.stdout",
                    },
                "file":
                    {
                        "class": "logging.FileHandler",
                        "formatter": "form02",
                        "level": file_loglevel,
                        "filename": log,
                        "mode": "w",  # Overwrite file if it exists
                    },
            },
        "loggers":
            {
                "":
                    {
                        "level": "NOTSET",
                        "handlers": ["console", "file"],
                    },
            },
        "disable_existing_loggers": False,
    }

    if extra_config is not None:
        log_config = update_dict(log_config, extra_config)

    logging.config.dictConfig(log_config)

def update_dict(d, u):
    """
    Recursively merge or update dict-like objects.
    i.e, change a value to a key that already exists or
    add a (key, value) that did not previously existed

    source: https://stackoverflow.com/questions/3232943/update-value-of-a-nested-dictionary-of-varying-depth

    :param dict d: Original dictionnary
    :param dict u: dictionnary of updates to apply to 'd'
    :return dict d: Return updated version of 'd'
    """

    for k, v in u.items():
        if isinstance(v, collections.Mapping):
            d[k] = update_dict(d.get(k, {}), v)
        else:
            d[k] = v
    return d

def analyse_apt_list(files):
    """
    Parse a list of APT files and store in a dictionnary the list of simulations corresponding to each APT file.

    The name of each .aptx file must be the ID of a STScI proposal. It must also correspond to a MIRI CAR.

    :param files: list of .aptx filenames
    :type files: list(str)

    :return: (APT_sim, mem_volume, time_volume) list of simulation parameter, prediction for memory and computation time for each car
    :rtype: (dict(car:list(dict(parameters))), dict(car:list(memory in GB)), dict(car:list(time in hours)))
    """
    cars = []
    fails = []

    APT_sim = {}
    for filename in files:
        basename = os.path.basename(filename)
        stsci_proposal_id = os.path.splitext(basename)[0]
        car_id = c.MIRI_CAR[stsci_proposal_id]

        simulations = parse_apt(filename)

        if simulations:
            cars.append(car_id)
            APT_sim[car_id] = simulations
        else:
            fails.append(car_id)

    if fails:
        LOG.error("The following files failed: {}".format(fails))

    (mem_volume, time_volume) = count_datavolume(APT_sim)

    return (APT_sim, mem_volume, time_volume)
