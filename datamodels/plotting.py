#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""

Plotting functions for the MIRI data model. Based on the visitor
design pattern.

:Reference:

The STScI jwst.datamodels documentation.
The miri.tools.miriplot module.

:History:

20 Jan 2011: Created for old data model
21 Jan 2013: Moved to new data model.
18 Jun 2013: Primary array is now obtained with the get_primary_array_name()
             method.
02 Jul 2013: Added the ability to get the X and Y axis labels from the
             data model.
05 Jul 2013: Check for data objects without a type.
30 Jul 2013: Check for attempts to plot null or scalar objects.
20 May 2014: Modified for jsonschema draft 4. Get data arrays, tables,
             title and unit strings, etc... using the data model methods
             rather than the schema structure.
11 Sep 2015: Issue a warning when an attempt is made to plot an empty
             data model. Close all plots when a plotting object is deleted.
04 Jan 2017: Make the module resilient to matplotlib import problems.
15 Jun 2017: Corrected issue with non-integer numpy indices (which fail
             under numpy 1.12.1).

@author: Steven Beard (UKATC)

"""
# This module is now converted to Python 3.


import warnings

# Import the miri.tools plotting module.
import miri.tools.miriplot as mplt

# List all classes and global functions here.
__all__ = ['DataModelVisitor', 'DataModelPlotVisitor']


class DataModelVisitor(object):
    """
    
    Class DataModelVisitor - defines a generic data product visitor
    class. A class of this type can be used to process a data product
    with a foreign algorithm without building a dependency into the
    data product itself. Its implementation is based on the visitor
    design pattern.
    
    :Interface:
    
    A visitor class is expected to implement the visit() method,
    which calls the visit_xxx method appropriate for the class
    of object provided in its parameter.
    
    """
    def visit(self, data_element, **kwargs):
        """
        
        Apply the algorithm known to a a visitor class to the given
        data element. This method switches control to the visit_xxx
        method appropriate for the class of object provided.
                 
        :Parameters:
    
        data_element: object
            The data element to which the algorithm is to be applied.
        \*\*kwargs
            Any unrecognised keywords are passed to the algorithm.
                
        :Returns:
        
        Whatever the visit_xxx method returns.
        
        """
        meth = None
        for cls in data_element.__class__.__mro__:
            meth_name = 'visit_' + cls.__name__
            meth = getattr(self, meth_name, None)
            if meth is not None:
                break
        if meth is None:
            meth = self.visit_generic
        return meth(data_element, **kwargs)
    
    def visit_generic(self, data_element):
        """
        
        The visit_xxx method called when the data element class is
        not recognised.
        
        """
        # Ignore unrecognised elements or raise an exception,
        # depending on how strictly you want to enforce the rules.
        raise TypeError("Unrecognised data structure element type: " + \
                        data_element.__class__.__name__)


class DataModelPlotVisitor(DataModelVisitor):
    """
    
    Class DataModelPlotting - defines a class which contains the
    knowledge of how to display a data structure. Based on the visitor
    design pattern.
    
    This particular class uses matplotlib. A different graphics
    package could be used by substituting this class with a different
    plotting class.
                             
    :Parameters:
    
    description: str, optional
        An optional string to be appended to the plot titles.
    
    :Interface:
    
    Contains a selection of visit_xxxx classes which implement the
    visitor design pattern interface.
    
    :Requires:
        
    miri.tools.miriplot
    matplotlib

    """
    def __init__(self, description=''):
        """
        
        Initialises the DataModelPlotVisitor class.
        
        Parameters: See class doc string.

        """
        # Is plotting available?
        self.available = mplt.test_available()
        
        # Initialise the plot status
        self.supertitle = ''
        self.description = description
        self.subrows = 1
        self.subcols = 1
        self.subno = 0
        self.nfigs = 1
        self.page = 0
        self.fig = None
        self.axis = None

    def __del__(self):
        """
        
        Tidies up the DataModelPlotVisitor class.
        
        If a plot page is still open it is closed.

        """
        # Close the plot
        if self.fig is not None:
            self._close_page()
        mplt.close()

    def _open_page(self):
        """
        
        Open a new plot page.
       
        """
        if self.subno > 0:
            self._close_page()
        self.page += 1
        self.subno = 0
        self.fig = mplt.new_figure(self.page, figsize=(10,10),
                                   stitle=self.supertitle)
        
    def _close_page(self, prompt=''):
        """
        
        Close the plot page.
        
        :Parameters:
        
        prompt: str, optional
            An optional string to prompt the user to close the plot 
            window.
       
        """
        if self.subno > 0:
            mplt.show_plot(prompt=prompt)
        self.fig = None
        self.axis = None

    def _new_subplot(self):
        """
        
        Create a new subplot.
       
        """
        if self.fig is None:
            self._open_page()
        elif self.subno >= self.subrows*self.subcols:
            self._open_page()
        self.subno += 1
        self.axis = mplt.add_subplot(self.fig, subrows=self.subrows,
                                     subcols=self.subcols, subno=self.subno)

    def visit_DataModel(self, data_element, **kwargs):
        """
        
        Apply the plotting algorithm to a any STScI DataModel object.
        This function interrogates the schema describing a STScI data
        object and plots the data arrays and data tables contained
        in the object.
        
        :Parameters:
        
        data_element: STScI DataModel
            The data structure to be plotted.
        \*\*kwargs
            Any unrecognised keywords are passed to the plotting function.
       
        """
        # The data structure must have a schema or the information
        # needed for the plot cannot be extracted.
        if not hasattr(data_element, "schema"):
            raise AttributeError("STScI DataModel object has no schema")
        
        if self.available:
            # Obtain a plot title from the schema model. Fall back on the
            # class name if the schema has no title. Add the description
            # supplied with the plot() function call.
            if 'title' in data_element.schema:
                self.supertitle = data_element.schema['title']
            else:
                self.supertitle = data_element.__class__.__name__
            if self.description:
                self.supertitle += " - %s" % self.description
    
            # Search the schema and compile a list of data arrays and
            # data tables belonging to this object.
            list_of_arrays = data_element.list_data_arrays()
            list_of_tables = data_element.list_data_tables()
            if not list_of_arrays and not list_of_tables:
                strg = "\nData model %s contains no plottable data." % \
                    data_element.__class__.__name__
                warnings.warn(strg)
                return
    
            # Subdivide the figures into a grid of subplots, one for each data
            # array and data table contained in the data product.
            # Obtain the dimensionality of the array either from the primary
            # array or from the first array in the list.
            # A 3-D data product needs each new plot on a different page.
            narrays = len(list_of_arrays)
            ntables = len(list_of_tables)
            nplots = narrays + ntables
            primary_name = data_element.get_primary_array_name()
            if hasattr(data_element, primary_name):
                primary_array = getattr(data_element, primary_name)
                firstdim = primary_array.ndim
            elif len(list_of_arrays) > 0:
                firstarray = list_of_arrays[0]
                first_array = getattr(data_element, firstarray)
                firstdim = first_array.ndim 
            else:
                firstdim = 0
            if firstdim < 3:
                (self.subrows, self.subcols, self.nfigs) = mplt.subdivide(nplots)
            else:
                self.subrows = 1
                self.subcols = 1
                self.nfigs = nplots
            
            # Plot the metadata
            self.visit_Meta(data_element.meta, **kwargs)
            
            # Plot each of the data arrays contained in the model.
            for array_name in list_of_arrays:
                if hasattr(data_element, 'dataname_to_hduname'):
                    name = data_element.dataname_to_hduname(array_name)
                    if name is None:
                        name = array_name
                else:
                    name = array_name
                if hasattr(data_element, 'get_data_units'):
                    units = data_element.get_data_units(array_name)
                else:
                    units = ""
                if hasattr(data_element, 'get_data_title'):
                    title = data_element.get_data_title(array_name)
                else:
                    title = ""
                xlabel = None
                ylabel = None
                if hasattr(data_element, 'get_data_axes'):
                    axes = data_element.get_data_axes(name)
                    if axes:
                        xlabel = axes[-1]
                        if len(axes) > 1:
                            ylabel = axes[-2]
    
                # If there is a "filled" version of the data array, plot
                # that in preference to the original array.
                filled_name = array_name + "_filled"
                if hasattr(data_element, filled_name):   
                    data_array = getattr(data_element, filled_name)
                    name += " (filled)"
                else:
                    data_array = getattr(data_element, array_name)
                self.visit_ndarray(data_array, name=name, units=units, title=title,
                                   xlabel=xlabel, ylabel=ylabel, **kwargs)
            
            # Plot each of the data tables contained in the model.
            for table_name in list_of_tables:
                if hasattr(data_element, 'dataname_to_hduname'):
                    name = data_element.dataname_to_hduname(table_name)
                    if name is None:
                        name = table_name
                else:
                    name = table_name
                if hasattr(data_element, 'get_data_title'):
                    title = data_element.get_data_title(table_name)
                else:
                    title = ""
                table = getattr(data_element, table_name)
                self.visit_recarray(table, name=name, title=title, **kwargs)
            self._close_page()
        else:
            strg = "MIRI plotting utility is not available "
            strg += "due to matplotlib import problem."
            warnings.warn(strg)

    def visit_Meta(self, data_element, **kwargs):
        """
        
        Apply the plotting algorithm to a Metadata object.
        
        :Parameters:
        
        data_element: Metadata
            The data element to be plotted.
        \*\*kwargs
            Any unrecognised keywords are passed to the plotting function.
       
        """       
        # For now, metadata is ignored.
        return
        
        # Metadata could be included using the code below (if space
        # for an extra subplot is made available).
#        if self.fig is None:
#            self._open_page()
#        self._new_subplot()
#        metastr = data_element.to_tree()
#        mplt.plot_text(metastr, title='Metadata', plotfig=self.fig,
#                       plotaxis=self.axis, **kwargs)

    def visit_ndarray(self, data_element, name="", units="", title="",
                      xlabel=None, ylabel=None, **kwargs):
        """
        
        Apply the plotting algorithm to a numpy ndarray object.
        This function plots a 1-D or 2-D data array within the next
        available subplot.
        A 3-D data array is plotted on a new page, with each plane
        plotted in a subplot.
        Data with dimensionality of 4 or above is plotted as a 3-D
        view.
        
        :Prior conditions:
        
        It is assumed this function is called in the same order in
        which the data arrays appear in the schema.
        
        :Parameters:
        
        data_element: numpy ndarray
            The data element to be plotted.
        name: str, optional
            The name of the data element, if any.
        units: str, optional
            The data units, if any.
        title: str, optional
            A title for the plot, if any.
        xlabel: str, optional
            A label for the X axis, if any. A default will be used if
            not specified.
        ylabel: str, optional
            A label for the Y axis, if any. A default will be used if
            not specified.
        \*\*kwargs
            Any unrecognised keywords are passed to the plotting function.
       
        """               
        if self.available:
            # Each plot is labelled with the name, type and unit of the data.
            tstrg = ""
            if title:
                tstrg += "%s: " % title
            if name:
                tstrg += "%s" % name
            if units:
                tstrg += " (%s)" % units
    
            # Null data cannot be plotted
            if data_element is None:
                strg = "\n***\'%s\' array is null and cannot be plotted." % tstrg
                warnings.warn(strg)
                return
     
            # 1-D and 2-D data can be plotted within a single subplot.
            # Higher dimensional data needs to start from a new page.
            oldshape = None
            if data_element.ndim <= 2:
                self._new_subplot()
            else:
                tstrg = self.supertitle + "\n" + tstrg
                oldsupertitle = self.supertitle
                self.supertitle = ''
                self._open_page()
                
                # Data with a dimensionality higher than 3-D is plotted
                # as a 3-D view.
                if data_element.ndim > 3:
                    olddim = data_element.ndim
                    oldshape = data_element.shape
                    slices = int(data_element.size / \
                        (oldshape[-2] * oldshape[-1]))
                    newshape = [slices, oldshape[-2], oldshape[-1]]
                    data_element.shape = newshape
                    tstrg += " (3-D view of %d-D data)" % olddim
                    
            if data_element.ndim < 1:
                # The object is a scalar
                mplt.plot_text(str(data_element), title=tstrg, plotfig=self.fig,
                               plotaxis=self.axis, **kwargs)
    
            elif len(data_element) < 1:
                # Empty data sets cannot be plotted
                mplt.plot_text('Empty array', title=tstrg, plotfig=self.fig,
                               plotaxis=self.axis, **kwargs)      
            elif data_element.ndim == 1:
                if not xlabel:
                    xlabel = 'Elements'
                mplt.plot_xy(None, data_element, xlabel=xlabel, ylabel=tstrg,
                             plotfig=self.fig, plotaxis=self.axis, **kwargs)
            else:
                mplt.plot_image(data_element, title=tstrg, xlabel=xlabel,
                                ylabel=ylabel, plotfig=self.fig,
                                plotaxis=self.axis, **kwargs)
                
            if data_element.ndim > 2:
                self.supertitle = oldsupertitle
                self._close_page()
                
            # Restore the previous shape of the data array.
            if oldshape is not None:
                data_element.shape = oldshape
        else:
            strg = "MIRI plotting utility is not available "
            strg += "due to matplotlib import problem."
            warnings.warn(strg)
 
    def visit_recarray(self, data_element, name="", title="", **kwargs):
        """
        
        Apply the plotting algorithm to a numpy recarray object.
        This function displays a data table as a textual list
        of columns.
        
        :Parameters:
        
        data_element: numpy recarray
            The data element to be plotted.
        name: str, optional
            The name of the data table, if any.
        title: str, optional
            A title for the plot, if any.
       
        """
        if self.available:
            # Plot the information contained in a data table as text.
            tstrg = ""
            if title:
                tstrg += "%s: " % title
            if name:
                tstrg += "%s" % name
    
            # Plot the contents of a data table as text
            self._new_subplot()
            tbstr = str(data_element)
            mplt.plot_text(tbstr, title=tstrg, plotfig=self.fig,
                           plotaxis=self.axis, **kwargs)
        else:
            strg = "MIRI plotting utility is not available "
            strg += "due to matplotlib import problem."
            warnings.warn(strg)


#
# A minimal test is run when this file is run as a main program.
# For a more substantial test see miri/datamodels/tests.
#
if __name__ == '__main__':
    print("Testing the plotting module.")
    # TBD
    print("Test finished.")
