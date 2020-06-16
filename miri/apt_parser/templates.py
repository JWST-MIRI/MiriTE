import re
import logging

from . import dithering
from . import constants as c

LOG = logging.getLogger('parse_apt.templates')


class Template(object):
    """
    Generic template object from which every class derive
    """

    def __init__(self, template, metadata=None):
        """

        :param template: MiriImaging observation XML object
        :type template: <class 'xml.etree.ElementTree.Element'>
        """
        self.template = template
        self.simulation_list = []
        self.metadata = metadata
        self.metadata["template"] = self.__class__.__name__
        self.parse_template()

    def extract_dither_id(self, text):
        """
        In observations, the dither is given like this: "Dither 1",
        this correspond to the 0-th element in the dither list.

        This function return 0 if given Dither 1, 1 if given Dither 2, etc...

        The function crash if it can't match the expected pattern

        :param str text: Text from filter_element.find("mi:Dither", parser.ns).text
        :return int: index of the corresponding dither in the dither list
        """

        regexp = "^Dither ([0-9]+)$"

        pattern_obj = re.compile(regexp)

        match = pattern_obj.match(text)

        if match:
            return int(match.group(1)) - 1
        else:
            raise ValueError("Regexp {} don't match in '{}'".format(regexp, text))

    def getobs(self):
        return self.simulation_list

    def parse_template(self):
        pass

    def _add_simulation(self, simulation):
        """
        Before adding the simulation to the list, will make a copy of common metadata from self.metadata

        :param simulation: Simulation dictionnary to add to the list of simulation for this Template
        :type simulation: dict()
        """

        simulation.update(self.metadata)
        self.simulation_list.append(simulation)

    def _extend_simulation_list(self, sim_list):
        """
        Similar to _add_simulation, but add several simulations at once

        :param sim_list: list of Simulation dictionnary to add to the list of simulation for this Template
        :type sim_list: list(dict())
        :return:
        :rtype:
        """

        for sim in sim_list:
            self._add_simulation(sim)


class MiriLRS(Template):
    NS = "mlrs"

    def parse_template(self):
        self.simulation_list = []

        subarray = self.template.find("{}:Subarray".format(self.NS), c.ns).text
        exposures = int(self.template.find("{}:Exposures".format(self.NS), c.ns).text)
        integrations = int(self.template.find("{}:Integrations".format(self.NS), c.ns).text)
        groups = int(self.template.find("{}:Groups".format(self.NS), c.ns).text)
        readout_pattern = self.template.find("{}:ReadoutPattern".format(self.NS), c.ns).text

        dither = self.template.find("{}:DitherType".format(self.NS), c.ns).text

        nb_dither_points = 1  # If no dithers
        if dither != 'NONE':
            try:
                nb_dither_points = c.DITHER_POINTS[dither]
            except KeyError:

                if dither == "MAPPING":
                    spatial_step = int(self.template.find("{}:NumberOfSpatialSteps".format(self.NS), c.ns).text)
                    spectral_step = int(self.template.find("{}:NumberOfSpectralSteps".format(self.NS), c.ns).text)
                    nb_dither_points = spatial_step * spectral_step
                else:
                    LOG.error("Unable to retrieve number of dither points from {}".format(dither))

        # Specific to my needs
        simulation = {
            "exposures": exposures,
            "ima_integrations": integrations,
            "ima_frames": groups,
            "subarray": subarray,
            "NDither": nb_dither_points,
            "readDetect": readout_pattern,
        }
        self._add_simulation(simulation)


class MiriImaging(Template):
    NS = "mi"

    def parse_template(self):
        self.simulation_list = []

        # Empty list if no dithers
        dithers = self.template.findall(".//{}:DitherSpecification".format(self.NS), c.ns)
        subarray = self.template.find("{}:Subarray".format(self.NS), c.ns).text

        filters = self.template.find("{}:Filters".format(self.NS), c.ns)

        for filter_element in filters:
            exposures = int(filter_element.find("{}:Exposures".format(self.NS), c.ns).text)
            integrations = int(filter_element.find("{}:Integrations".format(self.NS), c.ns).text)
            groups = int(filter_element.find("{}:Groups".format(self.NS), c.ns).text)
            filter_name = filter_element.find("{}:Filter".format(self.NS), c.ns).text
            readout_pattern = filter_element.find("{}:ReadoutPattern".format(self.NS), c.ns).text

            dither = filter_element.find("{}:Dither".format(self.NS), c.ns).text

            nb_dither_points = 1  # If no dithers
            if dither != 'None':
                dither_idx = self.extract_dither_id(dither)
                dither_el = dithers[dither_idx]

                nb_dither_points = dithering.parse_dither(dither_el, self.NS)

            # Specific to my needs
            simulation = {
                "exposures": exposures,
                "ima_integrations": integrations,
                "ima_frames": groups,
                "NDither": nb_dither_points,
                "subarray": subarray,
                "filter": filter_name,
                "readDetect": readout_pattern,
            }
            self._add_simulation(simulation)


class MiriExternalFlat(Template):
    NS = "mef"

    def parse_imager_exposure(self, exposure, common_metadata):
        """

        :param exposure: Exposure template
        :type exposure: mef:Exposure XML element
        :param dict common_metadata: metadata declared in parent element that need to be added
        :return: list of sub simulations
        :rtype:
        """

        exposures = int(exposure.find("{}:Exposures".format(self.NS), c.ns).text)
        integrations = int(exposure.find("{}:Integrations".format(self.NS), c.ns).text)
        groups = int(exposure.find("{}:Groups".format(self.NS), c.ns).text)
        filter_name = exposure.find("{}:Filter".format(self.NS), c.ns).text
        readout_pattern = exposure.find("{}:ReadoutPattern".format(self.NS), c.ns).text

        nb_dither = 1  # No dither by default
        if "NDither" in common_metadata:
            nb_dither = common_metadata["NDither"]

        # Specific to my needs
        simulation = {
            "exposures": exposures,
            "ima_integrations": integrations,
            "ima_frames": groups,
            "NDither": nb_dither,
            "filter": filter_name,
            "readDetect": readout_pattern,
        }

        simulation.update(common_metadata)

        return [simulation]

    def parse_mrs_exposure(self, exposure, common_metadata):
        """

        :param exposure: Exposure template
        :type exposure: mef:Exposure XML element
        :param dict common_metadata: metadata declared in parent element that need to be added
        :return: list of sub simulations
        :rtype:
        """

        exposures = int(exposure.find("{}:Exposures".format(self.NS), c.ns).text)
        readout_pattern = exposure.find("{}:ReadoutPatternLong".format(self.NS), c.ns).text

        disperser_long = exposure.find("{}:Wavelength1_4".format(self.NS), c.ns).text
        groups_long = int(exposure.find("{}:GroupsLong".format(self.NS), c.ns).text)
        integrations_long = int(exposure.find("{}:IntegrationsLong".format(self.NS), c.ns).text)

        disperser_short = exposure.find("{}:Wavelength2_3".format(self.NS), c.ns).text
        groups_short = int(exposure.find("{}:GroupsShort".format(self.NS), c.ns).text)
        integrations_short = int(exposure.find("{}:IntegrationsShort".format(self.NS), c.ns).text)

        # Since long and short can't be different I read only one
        # wavelength1_4 = exposure.find("{}:Wavelength1_4".format(self.NS), parser.ns).text
        # disperser = parser.DISPERSER[wavelength1_4]

        # Artificially add a subarray to predict memory and time consumption
        if disperser_long is not None and disperser_short is not None:
            detector = "BOTH"
            subarray = "BOTH"
        elif disperser_long is not None:
            detector = "LW"
            subarray = "FULL"
        elif disperser_short is not None:
            detector = "SW"
            subarray = "FULL"
        else:
            LOG.error("Can't determine detector (LW, SW or BOTH). disperser_long={}; disperser_short={}".format(
                disperser_long, disperser_short))

        if detector == "BOTH":
            if disperser_long != disperser_short:
                raise ValueError("[detector=BOTH] Incoherent values (disperser) for SW and LW")

            disperser = c.DISPERSER[disperser_long]
        elif detector == "LW":
            disperser = c.DISPERSER[disperser_long]
        elif detector == "SW":
            disperser = c.DISPERSER[disperser_short]

        nb_dither = 1  # No dither by default
        if "NDither" in common_metadata:
            nb_dither = common_metadata["NDither"]

        # Specific to my needs
        simulation = {
            "exposures": exposures,
            "LW_integrations": integrations_long,
            "LW_frames": groups_long,
            "SW_integrations": integrations_short,
            "SW_frames": groups_short,
            "subarray": subarray,
            "detector": detector,
            "disperser": disperser,
            "NDither": nb_dither,
            "readDetect": readout_pattern,
        }

        simulation.update(common_metadata)

        return [simulation]

    def parse_all_exposure(self, exposure, common_metadata):
        """

        :param exposure: Exposure template
        :type exposure: mef:Exposure XML element
        :param dict common_metadata: metadata declared in parent element that need to be added
        :return: list of sub simulations
        :rtype:
        """

        # IMA
        exposures = int(exposure.find("{}:Exposures".format(self.NS), c.ns).text)
        ima_integrations = int(exposure.find("{}:Integrations".format(self.NS), c.ns).text)
        ima_groups = int(exposure.find("{}:Groups".format(self.NS), c.ns).text)
        ima_filter = exposure.find("{}:Filter".format(self.NS), c.ns).text
        readout_pattern = exposure.find("{}:ReadoutPattern".format(self.NS), c.ns).text

        # Since long and short can't be different I read only one
        wavelength1_4 = exposure.find("{}:Wavelength1_4".format(self.NS), c.ns).text
        disperser = c.DISPERSER[wavelength1_4]

        readout_pattern_long = exposure.find("{}:ReadoutPatternLong".format(self.NS), c.ns).text
        groups_long = int(exposure.find("{}:GroupsLong".format(self.NS), c.ns).text)
        integrations_long = int(exposure.find("{}:IntegrationsLong".format(self.NS), c.ns).text)

        groups_short = int(exposure.find("{}:GroupsShort".format(self.NS), c.ns).text)
        integrations_short = int(exposure.find("{}:IntegrationsShort".format(self.NS), c.ns).text)

        nb_dither = 1  # No dither by default
        if "NDither" in common_metadata:
            nb_dither = common_metadata["NDither"]
        # Specific to my needs
        simulation = {
            "exposures": exposures,
            "filter": ima_filter,
            "ima_integrations": ima_integrations,
            "ima_frames": ima_groups,
            "LW_integrations": integrations_long,
            "LW_frames": groups_long,
            "SW_integrations": integrations_short,
            "SW_frames": groups_short,
            "disperser": disperser,
            "NDither": nb_dither,
            "readDetect": readout_pattern,
        }

        simulation.update(common_metadata)

        return [simulation]

    def parse_template(self):
        self.simulation_list = []

        detector = self.template.find("{}:Detector".format(self.NS), c.ns).text
        dither = self.template.find("{}:Dither".format(self.NS), c.ns).text

        exposure_list = self.template.find("{}:ExposureList".format(self.NS), c.ns)

        common_metadata = {
            "detector": detector,
            "dither_type": dither,
        }

        if dither == "true":
            dither_el = self.template.find("{}:DitherSpec".format(self.NS), c.ns)
            nb_dither_point = dithering.parse_dither(dither_el, self.NS)
            common_metadata["NDither"] = nb_dither_point

        if detector == "IMAGER":
            self.parse_exposure = self.parse_imager_exposure
            subarray = self.template.find("{}:Subarray".format(self.NS), c.ns).text
            common_metadata["subarray"] = subarray
        elif detector == "MRS":
            self.parse_exposure = self.parse_mrs_exposure
        elif detector == "ALL":
            # Dual mode MRS+IMA
            subarray = self.template.find("{}:Subarray".format(self.NS), c.ns).text
            common_metadata["subarray"] = subarray

            self.parse_exposure = self.parse_all_exposure
        else:
            raise ValueError("Unknown detector value ({})".format(detector))

        for exp_el in exposure_list:
            sim_list = self.parse_exposure(exp_el, common_metadata)
            self._extend_simulation_list(sim_list)


class MiriDark(Template):
    NS = "md"

    def parse_imager_exposure(self, exposure, common_metadata):
        """

        :param exposure: Exposure template
        :type exposure: mef:Exposure XML element
        :param dict common_metadata: metadata declared in parent element that need to be added
        :return: list of sub simulations
        :rtype:
        """

        exposures = int(exposure.find("{}:Exposures".format(self.NS), c.ns).text)
        integrations = int(exposure.find("{}:Integrations".format(self.NS), c.ns).text)
        groups = int(exposure.find("{}:Groups".format(self.NS), c.ns).text)
        readout_pattern = exposure.find("{}:ReadoutPattern".format(self.NS), c.ns).text

        nb_dither = 1  # No dither by default

        # Specific to my needs
        simulation = {
            "exposures": exposures,
            "ima_integrations": integrations,
            "ima_frames": groups,
            "NDither": nb_dither,
            "readDetect": readout_pattern,
        }

        simulation.update(common_metadata)

        return [simulation]

    def parse_mrs_exposure(self, exposure, common_metadata):
        """

        :param exposure: Exposure template
        :type exposure: mef:Exposure XML element
        :param dict common_metadata: metadata declared in parent element that need to be added
        :return: list of sub simulations
        :rtype:
        """

        exposures = int(exposure.find("{}:Exposures".format(self.NS), c.ns).text)
        readout_pattern = exposure.find("{}:ReadoutPattern".format(self.NS), c.ns).text

        groups_long = int(exposure.find("{}:GroupsLong".format(self.NS), c.ns).text)
        integrations_long = int(exposure.find("{}:IntegrationsLong".format(self.NS), c.ns).text)

        groups_short = int(exposure.find("{}:GroupsShort".format(self.NS), c.ns).text)
        integrations_short = int(exposure.find("{}:IntegrationsShort".format(self.NS), c.ns).text)

        # Artificially add a subarray to predict memory and time consumption
        if groups_long is not None and groups_short is not None:
            detector = "BOTH"
            subarray = "BOTH"
        elif groups_long is not None:
            detector = "LW"
            subarray = "FULL"
        elif groups_short is not None:
            detector = "SW"
            subarray = "FULL"
        else:
            LOG.error("Can't determine detector (LW, SW or BOTH).")

        nb_dither = 1  # No dither by default

        # Specific to my needs
        simulation = {
            "exposures": exposures,
            "LW_integrations": integrations_long,
            "LW_frames": groups_long,
            "SW_integrations": integrations_short,
            "SW_frames": groups_short,
            "subarray": subarray,
            "detector": detector,
            "NDither": nb_dither,
            "readDetect": readout_pattern,
        }

        simulation.update(common_metadata)

        return [simulation]

    def parse_all_exposure(self, exposure, common_metadata):
        """

        :param exposure: Exposure template
        :type exposure: mef:Exposure XML element
        :param dict common_metadata: metadata declared in parent element that need to be added
        :return: list of sub simulations
        :rtype:
        """

        # IMA
        exposures = int(exposure.find("{}:Exposures".format(self.NS), c.ns).text)
        ima_integrations = int(exposure.find("{}:Integrations".format(self.NS), c.ns).text)
        ima_groups = int(exposure.find("{}:Groups".format(self.NS), c.ns).text)
        readout_pattern = exposure.find("{}:ReadoutPattern".format(self.NS), c.ns).text

        readout_pattern_long = exposure.find("{}:ReadoutPatternLong".format(self.NS), c.ns).text
        groups_long = int(exposure.find("{}:GroupsLong".format(self.NS), c.ns).text)
        integrations_long = int(exposure.find("{}:IntegrationsLong".format(self.NS), c.ns).text)

        readout_pattern_short = exposure.find("{}:ReadoutPatternShort".format(self.NS), c.ns).text
        groups_short = int(exposure.find("{}:GroupsShort".format(self.NS), c.ns).text)
        integrations_short = int(exposure.find("{}:IntegrationsShort".format(self.NS), c.ns).text)

        nb_dither = 1  # No dither by default
        if "NDither" in common_metadata:
            nb_dither = common_metadata["NDither"]
        # Specific to my needs
        simulation = {
            "exposures": exposures,
            "ima_integrations": ima_integrations,
            "ima_frames": ima_groups,
            "LW_integrations": integrations_long,
            "LW_frames": groups_long,
            "SW_integrations": integrations_short,
            "SW_frames": groups_short,
            "NDither": nb_dither,
            "readDetect": readout_pattern,
            "LW_readDetect": readout_pattern_long,
            "SW_readDetect": readout_pattern_short,
        }

        simulation.update(common_metadata)

        return [simulation]

    def parse_template(self):
        self.simulation_list = []

        detector = self.template.find("{}:Detector".format(self.NS), c.ns).text

        exposure_list = self.template.find("{}:Filters".format(self.NS), c.ns)

        common_metadata = {
            "detector": detector,
        }

        if detector == "IMAGER":
            self.parse_exposure = self.parse_imager_exposure
            subarray = self.template.find("{}:Subarray".format(self.NS), c.ns).text
            common_metadata["subarray"] = subarray
        elif detector == "MRS":
            self.parse_exposure = self.parse_mrs_exposure
        elif detector == "ALL":
            # Dual mode MRS+IMA
            subarray = self.template.find("{}:Subarray".format(self.NS), c.ns).text
            common_metadata["subarray"] = subarray

            self.parse_exposure = self.parse_all_exposure
        else:
            raise ValueError("Unknown detector value ({})".format(detector))

        for exp_el in exposure_list:
            sim_list = self.parse_exposure(exp_el, common_metadata)
            self._extend_simulation_list(sim_list)


class MiriMRS(Template):
    NS = "mmrs"

    def parse_mrs_exposure(self, exposure, common_metadata, dithers=None):
        """

        :param exposure: Exposure template
        :type exposure: mef:Exposure XML element
        :param dict common_metadata: metadata declared in parent element that need to be added
        :param dithers: [optional] Dithers template, contain a list of dither_type object
        :type dithers: mmrs:Dithers XML element
        :return: list of sub simulations
        :rtype:
        """

        exposures = int(exposure.find("{}:Exposures".format(self.NS), c.ns).text)
        readout_pattern = exposure.find("{}:ReadoutPatternLong".format(self.NS), c.ns).text

        disperser_band = exposure.find("{}:Wavelength".format(self.NS), c.ns).text
        disperser = c.DISPERSER[disperser_band]

        groups_long = int(exposure.find("{}:GroupsLong".format(self.NS), c.ns).text)
        integrations_long = int(exposure.find("{}:IntegrationsLong".format(self.NS), c.ns).text)

        groups_short = int(exposure.find("{}:GroupsShort".format(self.NS), c.ns).text)
        integrations_short = int(exposure.find("{}:IntegrationsShort".format(self.NS), c.ns).text)

        # Artificially add a subarray to predict memory and time consumption
        if groups_long is not None and groups_short is not None:
            detector = "BOTH"
            subarray = "BOTH"
        elif groups_long is not None:
            detector = "LW"
            subarray = "FULL"
        elif groups_short is not None:
            detector = "SW"
            subarray = "FULL"
        else:
            LOG.error("Can't determine detector (LW, SW or BOTH).")

        dither = exposure.find("{}:Dither".format(self.NS), c.ns).text

        nb_dither_points = 1  # No dither by default
        if dither != 'None':
            dither_idx = self.extract_dither_id(dither)
            dither_el = dithers[dither_idx]
            nb_dither_points = dithering.parse_dither(dither_el, self.NS)

        # Specific to my needs
        simulation = {
            "exposures": exposures,
            "LW_integrations": integrations_long,
            "LW_frames": groups_long,
            "SW_integrations": integrations_short,
            "SW_frames": groups_short,
            "subarray": subarray,
            "detector": detector,
            "disperser": disperser,
            "NDither": nb_dither_points,
            "readDetect": readout_pattern,
        }

        simulation.update(common_metadata)

        return [simulation]

    def parse_all_exposure(self, exposure, common_metadata, dithers=None):
        """

        :param exposure: Exposure template
        :type exposure: mef:Exposure XML element
        :param dict common_metadata: metadata declared in parent element that need to be added
        :param dithers: [optional] Dithers template, contain a list of dither_type object
        :type dithers: mmrs:Dithers XML element
        :return: list of sub simulations
        :rtype:
        """

        # IMA
        exposures = int(exposure.find("{}:Exposures".format(self.NS), c.ns).text)
        ima_integrations = int(exposure.find("{}:Integrations".format(self.NS), c.ns).text)
        ima_groups = int(exposure.find("{}:Groups".format(self.NS), c.ns).text)
        ima_filter = exposure.find("{}:Filter".format(self.NS), c.ns).text
        readout_pattern = exposure.find("{}:ReadoutPattern".format(self.NS), c.ns).text

        # Since long and short can't be different I read only one
        wavelength = exposure.find("{}:Wavelength".format(self.NS), c.ns).text
        disperser = c.DISPERSER[wavelength]

        readout_pattern_long = exposure.find("{}:ReadoutPatternLong".format(self.NS), c.ns).text
        groups_long = int(exposure.find("{}:GroupsLong".format(self.NS), c.ns).text)
        integrations_long = int(exposure.find("{}:IntegrationsLong".format(self.NS), c.ns).text)

        groups_short = int(exposure.find("{}:GroupsShort".format(self.NS), c.ns).text)
        integrations_short = int(exposure.find("{}:IntegrationsShort".format(self.NS), c.ns).text)

        dither = exposure.find("{}:Dither".format(self.NS), c.ns).text

        nb_dither_points = 1  # No dither by default
        if dither != 'None':
            dither_idx = self.extract_dither_id(dither)
            dither_el = dithers[dither_idx]
            nb_dither_points = dithering.parse_dither(dither_el, self.NS)

        # Specific to my needs
        simulation = {
            "exposures": exposures,
            "filter": ima_filter,
            "ima_integrations": ima_integrations,
            "ima_frames": ima_groups,
            "LW_integrations": integrations_long,
            "LW_frames": groups_long,
            "SW_integrations": integrations_short,
            "SW_frames": groups_short,
            "disperser": disperser,
            "NDither": nb_dither_points,
            "readDetect": readout_pattern,
        }

        simulation.update(common_metadata)

        return [simulation]

    def parse_template(self):
        self.simulation_list = []

        detector = self.template.find("{}:Detector".format(self.NS), c.ns).text
        dithers = self.template.find("{}:Dithers".format(self.NS), c.ns)
        simultaneous_imaging = self.template.find("{}:SimultaneousImaging".format(self.NS), c.ns).text

        exposure_list = self.template.find("{}:ExposureList".format(self.NS), c.ns)

        common_metadata = {
            "detector": detector,
            "simultaneous_imaging": simultaneous_imaging,
        }

        if detector == "MRS":
            self.parse_exposure = self.parse_mrs_exposure
        elif detector == "ALL":
            # Dual mode MRS+IMA
            subarray = self.template.find("{}:Subarray".format(self.NS), c.ns).text
            primary_channel = self.template.find("{}:PrimaryChannel".format(self.NS), c.ns).text
            common_metadata["subarray"] = subarray
            common_metadata["primary_channel"] = primary_channel

            self.parse_exposure = self.parse_all_exposure
        else:
            raise ValueError("Unknown detector value ({})".format(detector))

        for exp_el in exposure_list:
            sim_list = self.parse_exposure(exp_el, common_metadata, dithers=dithers)
            self._extend_simulation_list(sim_list)


class MiriCoron(Template):
    NS = "mc"

    def parse_template(self):
        self.simulation_list = []

        # Empty list if no dithers
        dither = self.template.find("{}:Dither".format(self.NS), c.ns).text

        nb_dither_points = 1  # If no dithers
        if dither != 'NONE':
            try:
                nb_dither_points = c.DITHER_POINTS[dither]
            except KeyError:
                LOG.error("Unknown Dither type: {}".format(dither))

        filters = self.template.find("{}:Filters".format(self.NS), c.ns)

        for filter_element in filters:
            exposures = int(filter_element.find("{}:Exposures".format(self.NS), c.ns).text)
            integrations = int(filter_element.find("{}:Integrations".format(self.NS), c.ns).text)
            groups = int(filter_element.find("{}:Groups".format(self.NS), c.ns).text)
            filter_name = filter_element.find("{}:Filter".format(self.NS), c.ns).text
            mask = filter_element.find("{}:Mask".format(self.NS), c.ns).text
            readout_pattern = filter_element.find("{}:ReadoutPattern".format(self.NS), c.ns).text

            # Specific to my needs
            simulation = {
                "exposures": exposures,
                "ima_integrations": integrations,
                "ima_frames": groups,
                "NDither": nb_dither_points,
                "subarray": mask,
                "filter": filter_name,
                "readDetect": readout_pattern,
            }
            self._add_simulation(simulation)


class MiriCpc(Template):
    NS = "mcpc"

    def parse_template(self):
        self.simulation_list = []

        nb_dither_points = 1  # no dither

        subarray = self.template.find("{}:Subarray".format(self.NS), c.ns).text
        filters = self.template.find("{}:Filters".format(self.NS), c.ns)

        for filter_element in filters:
            exposures = int(filter_element.find("{}:Exposures".format(self.NS), c.ns).text)
            integrations = int(filter_element.find("{}:Integrations".format(self.NS), c.ns).text)
            groups = int(filter_element.find("{}:Groups".format(self.NS), c.ns).text)
            filter_name = filter_element.find("{}:Filter".format(self.NS), c.ns).text
            readout_pattern = filter_element.find("{}:ReadoutPattern".format(self.NS), c.ns).text

            # Specific to my needs
            simulation = {
                "exposures": exposures,
                "ima_integrations": integrations,
                "ima_frames": groups,
                "NDither": nb_dither_points,
                "subarray": subarray,
                "filter": filter_name,
                "readDetect": readout_pattern,
            }
            self._add_simulation(simulation)


class MiriMRSCrossGratingEngineering(Template):
    NS = "mmrscge"

    def parse_template(self):
        self.simulation_list = []

        exposure_list = self.template.find("{}:ExposureList".format(self.NS), c.ns)

        for exp_el in exposure_list:
            exposures = int(exp_el.find("{}:Exposures".format(self.NS), c.ns).text)

            groups_long = int(exp_el.find("{}:GroupsLong".format(self.NS), c.ns).text)
            integrations_long = int(exp_el.find("{}:IntegrationsLong".format(self.NS), c.ns).text)

            groups_short = int(exp_el.find("{}:GroupsShort".format(self.NS), c.ns).text)
            integrations_short = int(exp_el.find("{}:IntegrationsShort".format(self.NS), c.ns).text)

            # Specific to my needs
            simulation = {
                "exposures": exposures,
                "LW_integrations": integrations_long,
                "LW_frames": groups_long,
                "SW_integrations": integrations_short,
                "SW_frames": groups_short,
                "subarray": "FULL",
                "detector": "BOTH",
                "NDither": 1,
            }

            self._add_simulation(simulation)


templates = {"MiriImaging": MiriImaging, "MiriLRS": MiriLRS, "MiriExternalFlat": MiriExternalFlat, "MiriMRS": MiriMRS,
             "MiriDark": MiriDark, "MiriCoron": MiriCoron,
             "MiriCpc": MiriCpc, "MiriMRSCrossGratingEngineering": MiriMRSCrossGratingEngineering}
