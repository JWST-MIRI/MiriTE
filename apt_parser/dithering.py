import logging
LOG = logging.getLogger('parse_apt.dithering')

from . import constants as c


def parse_dither(dither_el, ns):
    """

    :param dither_el: mef:DitherSpec
    :type dither_el:
    :param ns: key for Namespace of the XML object (dictionary is in constants.py)
    :type ns: key
    :return:
    :rtype:
    """

    # For now we only want the total number of dithers
    # MIRISim only have 1 dither pattern implemented, not the others
    el_points = dither_el.find("./{}:Points".format(ns), c.ns)
    el_nb_points = dither_el.find("./{}:NumberOfPoints".format(ns), c.ns)
    dither_type = dither_el.find("./{}:DitherType".format(ns), c.ns).text

    if el_points is not None:
        point_detail = el_points.text
        if "," in point_detail:
            tmp = point_detail.split(",")
            nb_dither_points = 0
            for el in tmp:
                if "-" in el:
                    start, stop = el.split("-")
                    start = int(start)
                    stop = int(stop)
                    nb_dither_points += stop - start + 1
                else:
                    nb_dither_points += 1
        else:
            nb_dither_points = int(el_points.text)
    elif el_nb_points is not None:
        nb_dither_points = int(el_nb_points.text)
    else:
        try:
            nb_dither_points = c.DITHER_POINTS[dither_type]
        except KeyError:
            LOG.error("Unable to retrieve number of dither points from {}".format(dither_type))

    nb_sets = dither_el.find("./{}:NumberOfSets".format(ns), c.ns)
    if nb_sets is not None:
        nb_dither_points *= int(nb_sets.text)

    return nb_dither_points