import os.path
import asdf.util
import asdf.constants
from asdf.tags.core import HistoryEntry
from asdf.extension import AsdfExtension
from  asdf import AsdfFile
from asdf import schema as asdf_schema
from jwst import datamodels


SCHEMA_PATH = os.path.abspath(
    os.path.join(os.path.dirname(__file__), 'schemas'))


URL_PREFIX = asdf.constants.STSCI_SCHEMA_URI_BASE + 'miri_datamodel/'


METASCHEMA_PATH = os.path.abspath(
    os.path.join(os.path.dirname(datamodels.__file__), 'metaschema'))


class MIRIExtension(AsdfExtension):
    """                                                                                                      
                                                                                                             
    A class which defines an extension to the ASDF URI mapping                                               
    which recognises the URI 'http://stsci.edu/schemas/miri_datamodel'                                       
    as a reference to the schemas belonging to the miri.datamodels                                           
    package.                                                                                                 
                                                                                                             
    """

    @property
    def types(self):
        return []

    @property
    def tag_mapping(self):
        return []

    @property
    def url_mapping(self):
        # Define the uri and directory containing the additional schemas                                     
        return [
            (URL_PREFIX, asdf.util.filepath_to_url(SCHEMA_PATH) + '/{url_suffix}.yaml'),
            ('http://stsci.edu/schemas/fits-schema/',
             asdf.util.filepath_to_url(METASCHEMA_PATH) + '/{url_suffix}.yaml'),
        ]

