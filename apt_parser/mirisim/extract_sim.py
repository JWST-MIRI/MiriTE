"""
From the output of apt_parser, prepare dictionnaries for MIRISim. Each dictionnary represent one simulation
"""

from . import constants

import logging
import os

LOG = logging.getLogger('parse_apt.mirisim.extract_sim')

__all__ = ["get_imsim_parameters", "get_mrs_parameters", "get_lrs_parameters", "get_simulations"]

def get_imsim_parameters(observation):
    """
    Take an APT observation dictionnary and construct the MIRISim dictionnary of parameters needed to run a
    MIRI Imager simulation using MIRISim

    Note that for a parallel observation, one needs to run both function on the same observation to retrieve all
    needed simulations.

    :param observation: List of observations
    :return: dict of MIRISim simulations parameters (use SimConfig.makeSim(**dict) to create the object)
    """
    sim = constants.MIRISIM_DEFAULT.copy()
    sim["name"] = sim["name"] = "obs{}_IMA_{}".format(observation["obs_id"], observation["filter"])
    sim["POP"] = 'IMA'
    if "simultaneous_imaging" in observation:
        sim["DitherPat"] = "mrs_recommended_dither.dat"
    else:
        sim["DitherPat"] = "ima_recommended_dither.dat"

    sim["ConfigPath"] = "IMA_{}".format(observation["subarray"])
    sim["ima_exposures"] = observation["exposures"]

    if observation["NDither"] > 1:
        sim["Dither"] = True
    else:
        sim["Dither"] = False


    sim["readDetect"] = observation["subarray"]
    sim["ima_mode"] = observation["readDetect"]

    keys = ["filter", "NDither", "ima_frames", "ima_integrations"]
    for key in keys:
        sim[key] = observation[key]

    return sim


def get_lrs_parameters(observation):
    """
    Take an APT observation dictionnary and construct the MIRISim dictionnary of parameters needed to run a
    MIRI Imager simulation using MIRISim

    Note that for a parallel observation, one needs to run both function on the same observation to retrieve all
    needed simulations.

    :param observation: List of observations
    :return: dict of MIRISim simulations parameters (use SimConfig.makeSim(**dict) to create the object)
    """
    sim = constants.MIRISIM_DEFAULT.copy()
    sim["POP"] = 'IMA'
    sim["DitherPat"] = "lrs_recommended_dither.dat"
    sim["filter"] = "P750L"

    if observation["subarray"] == "FULL":
        sim["ConfigPath"] = "LRS_SLIT"
    else:
        sim["ConfigPath"] = "LRS_SLITLESS"

    sim["name"] = sim["name"] = "obs{}_{}".format(observation["obs_id"], sim["ConfigPath"])

    if observation["NDither"] > 1:
        sim["Dither"] = True
    else:
        sim["Dither"] = False


    sim["readDetect"] = observation["subarray"]
    sim["ima_mode"] = observation["readDetect"]
    sim["ima_exposures"] = observation["exposures"]

    keys = ["NDither", "ima_frames", "ima_integrations"]
    for key in keys:
        sim[key] = observation[key]

    return sim


def get_mrs_parameters(observation):
    """
    Take an APT observation dictionnary and construct the MIRISim dictionnary of parameters needed to run a
    MIRI MRS simulation using MIRISim

    Note that for a parallel observation, one needs to run both function on the same observation to retrieve all
    needed simulations.

    :param observation: List of observations
    :return: dict of MIRISim simulations parameters (use SimConfig.makeSim(**dict) to create the object)
    """

    sim = constants.MIRISIM_DEFAULT.copy()
    sim["name"] = "obs{}_MRS_{}".format(observation["obs_id"], observation["disperser"])
    sim["POP"] = 'MRS'
    sim["DitherPat"] = "mrs_recommended_dither.dat"

    if observation["detector"] in ["BOTH", "ALL", "MRS"]:
        sim["detector"] = "BOTH"
    else:
        sim["detector"] = observation["detector"]

    # MIRISim at the time do not take different values for the 2 MRS detectors.
    # Warning, we force to use the Long wavelength values.
    sim["mrs_exposures"] = observation["exposures"]
    sim["mrs_integrations"] = observation["LW_integrations"]
    sim["mrs_frames"] = observation["LW_frames"]

    if observation["NDither"] > 1:
        sim["Dither"] = True
    else:
        sim["Dither"] = False

    try:
        primary_channel = constants.PRIMARY_CHANNEL[observation["primary_channel"]]
    except KeyError:
        # Template MiriMRS doesn't have a primary_channel keyword.
        primary_channel = 1

    sim["ConfigPath"] = "MRS_{}{}".format(primary_channel, observation["disperser"])
    sim["mrs_mode"] = observation["readDetect"]

    keys = ["NDither", "disperser"]
    for key in keys:
        sim[key] = observation[key]


    return sim


def get_simulations(observations):
    """
    Take the output of the APT parser and duplicate simulations when needed to make individual separated simulations

    For instance, a parallel observation for MRS will generate an imager simulation, and as many as necessary MRS simulation
    for the required spectrum coverage.


    :param observations:
    :return: simulations
    """

    simulations = []
    for obs in observations:
        template = obs["template"]

        # Skipping templates impossible to do in MIRISim
        if template in constants.NON_SUPPORTED_TEMPLATES:
            LOG.warning("{} not supported by MIRISim. Skipping {} obs {}".format(template, os.path.basename(obs["filename"]), obs["obs_id"]))
            continue

        if not "detector" in obs or obs["simultaneous_imaging"] == "YES":
            if template == "MiriLRS":
                sim = get_lrs_parameters(obs)
            else:
                filter = obs["filter"]
                if filter == "FND":
                    LOG.warning("{} not supported by MIRISim. Skipping {} obs {}".format(filter, os.path.basename(
                        obs["filename"]), obs["obs_id"]))
                    continue

                sim = get_imsim_parameters(obs)
            simulations.append(sim)

        if "detector" in obs:
            # Multiple observations needs to be done
            sim = get_mrs_parameters(obs)
            simulations.append(sim)

    return simulations
