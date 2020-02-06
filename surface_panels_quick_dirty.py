"""quick workaround for comparing variables from different WRF runs

To do this right (TM) I need to figure out if/how I can generate
panel.depends decorators dynamically.  I haven't figured out the
syntax for that yet and I need these plots for the group meeting this
afternoon, so I'm going with copy/paste for the moment.

"""

import panel as pn
import geoviews_tools as gt


def three_panel_quadmesh_compare_vertical_var(ds, varname, cmap='RdBu'):
    """three-panel WRF variable comparison with sliders for z, time

    Create a three-panel plot showing values for W vertical wind
    velocity (W) with sliders to select vertical level and time stamp.

    Display a contour plot to show height above ground of the
    currently displayed vertical level.

    The three panels show values for the control run, Yatir run, and
    control - yatir difference.

    """

    hour_select = pn.widgets.IntSlider(start=0, end=24, value=9, name='Hour',
                                       orientation='vertical',
                                       direction='rtl')
    var_varies_vertically = len(gt.get_vdim(ds, varname)) > 0

    vdim = gt.get_vdim(ds, varname)

    if 'stag' in vdim:
        agl_var = 'height_agl_stag'
        zmax = ds[vdim].size
    elif vdim == []:   # this variable does not vary vertically
        zmax = 0
    else:
        agl_var = 'height_agl'
        zmax = ds[vdim].size

    z_select = pn.widgets.IntSlider(start=0, end=zmax, value=1,
                                    name='vertical level',
                                    orientation='vertical',
                                    direction='rtl',
                                    disabled=(var_varies_vertically is False))

    # bounds for the figures in fraction of the panel,
    fig_bounds = (0.2, 0.2, 0.8, 0.8)

    @pn.depends(hour_select, z_select)
    def get_contour_agl(hour_select, z_select):
        """create contours of height above ground level
        """
        # zstag: height of staggered Z levels, calucated by wrf-python
        # ter: height of terrain (meters above sea level)
        # calculate staggered Z level height above ground level (agl)
        agl_contour = ds[agl_var].sel({'WRFrun': 'control',
                                       'hour': hour_select,
                                       vdim: z_select}). hvplot.contour(
                                           x='XLONG',
                                           y='XLAT',
                                           z=agl_var,
                                           title='WRF height AGL').opts(
                                               fig_bounds=fig_bounds)
        return(agl_contour)

    @pn.depends(hour_select, z_select)
    def get_quadmesh_control(hour_select, z_select):
        """
        """
        # if var_varies_vertically:
        #     idx = {'WRFrun': 'control',
        #            'hour': hour_select,
        #            vdim: z_select}
        # else:
        #     idx = {'WRFrun': 'control',
        #            'hour': hour_select}
        vdim = gt.get_vdim(ds, varname)
        vmin, vmax = gt.get_min_max(ds, varname, hour_select, z_select)
        qm = ds[varname].sel({'WRFrun': 'control',
                              'hour': hour_select,
                              vdim: z_select}).hvplot.quadmesh(
                                  x='XLONG',
                                  y='XLAT',
                                  z=varname,
                                  title='control',
                                  clim=(vmin, vmax),
                                  cmap=cmap).opts(
                                      fig_bounds=fig_bounds)
        return(qm)

    @pn.depends(hour_select, z_select)
    def get_quadmesh_yatir(hour_select, z_select):
        """
        """
        vdim = gt.get_vdim(ds, varname)
        vmin, vmax = gt.get_min_max(ds, varname, hour_select, z_select)
        qm = ds[varname].sel({'WRFrun': 'yatir',
                              'hour': hour_select,
                              vdim: z_select}).hvplot.quadmesh(
                                  x='XLONG',
                                  y='XLAT',
                                  z=varname,
                                  title='Yatir',
                                  cmap=cmap,
                                  clim=(vmin, vmax)).opts(
                                      fig_bounds=fig_bounds)
        return(qm)

    @pn.depends(hour_select, z_select)
    def get_quadmesh_diff(hour_select, z_select):
        """
        """
        vdim = gt.get_vdim(ds, varname)
        vmin, vmax = gt.get_min_max(ds, varname, hour_select, z_select)
        qm = ds[varname].sel({'WRFrun': 'control - Yatir',
                              'hour': hour_select,
                              vdim: z_select}).hvplot.quadmesh(
                                  x='XLONG',
                                  y='XLAT',
                                  z=varname,
                                  #clim=(vmin, vmax),
                                  symmetric=True,
                                  cmap='RdBu',
                                  title='control - Yatir').opts(
                                      fig_bounds=fig_bounds)
        return(qm)

    main_title = '## ' + ds[varname].long_name
    the_plot = pn.Column(pn.Row(pn.pane.Markdown(main_title)),
                         pn.Row(get_quadmesh_control, get_quadmesh_yatir,
                                hour_select, z_select),
                         pn.Row(get_quadmesh_diff, get_contour_agl))
    return(the_plot)


def three_panel_quadmesh_compare_surface_var(ds, varname, cmap='RdBu'):
    """three-panel WRF variable comparison with slider for time

    Create a three-panel plot showing values for a WRF variable) with
    one slider to select time stamp.
    The three panels show values for the control run, Yatir run, and
    control - yatir difference.

    """

    hour_select = pn.widgets.IntSlider(start=0, end=24, value=9, name='Hour',
                                       orientation='vertical',
                                       direction='rtl')
    z_select = None  # this is a surface variable
    # bounds for the figures in fraction of the panel,
    fig_bounds = (0.2, 0.2, 0.8, 0.8)

    @pn.depends(hour_select)
    def get_quadmesh_control(hour_select):
        """
        """
        vmin, vmax = gt.get_min_max(ds, varname, hour_select, z_select)
        qm = ds[varname].sel({'WRFrun': 'control',
                              'hour': hour_select}).hvplot.quadmesh(
                                  x='XLONG',
                                  y='XLAT',
                                  z=varname,
                                  title='control',
                                  clim=(vmin, vmax),
                                  cmap=cmap).opts(
                                      fig_bounds=fig_bounds)
        return(qm)

    @pn.depends(hour_select)
    def get_quadmesh_yatir(hour_select):
        """
        """
        vmin, vmax = gt.get_min_max(ds, varname, hour_select, z_select)
        qm = ds[varname].sel({'WRFrun': 'yatir',
                              'hour': hour_select}).hvplot.quadmesh(
                                  x='XLONG',
                                  y='XLAT',
                                  z=varname,
                                  title='Yatir',
                                  cmap=cmap,
                                  clim=(vmin, vmax)).opts(
                                      fig_bounds=fig_bounds)
        return(qm)

    @pn.depends(hour_select)
    def get_quadmesh_diff(hour_select):
        """
        """
        vmin, vmax = gt.get_min_max(ds, varname, hour_select, z_select)
        qm = ds[varname].sel({'WRFrun': 'control - Yatir',
                              'hour': hour_select}).hvplot.quadmesh(
                                  x='XLONG',
                                  y='XLAT',
                                  z=varname,
                                  #clim=(vmin, vmax),
                                  symmetric=True,
                                  cmap='RdBu',
                                  title='control - Yatir').opts(
                                      fig_bounds=fig_bounds)
        return(qm)

    main_title = '## ' + ds[varname].long_name
    the_plot = pn.Column(pn.Row(pn.pane.Markdown(main_title)),
                         pn.Row(get_quadmesh_control, get_quadmesh_yatir,
                                hour_select),
                         pn.Row(get_quadmesh_diff))
    return(the_plot)
