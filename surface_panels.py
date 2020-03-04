import panel as pn
import numpy as np


def quadmesh_compare_vertical_var(ds,
                                  varname,
                                  cmap='RdBu',
                                  vdim='num_metgrid_levels'):
    """WRF variable comparison with sliders for z, time

    Create a plot showing values for a WRF variable with sliders to
    select vertical level and time stamp.

    """

    hour_select = pn.widgets.IntSlider(start=0,
                                       end=1, #ds['Time'].values.max(),
                                       value=0,
                                       name='Time',
                                       orientation='vertical',
                                       direction='rtl')
    if vdim is None:
        zmax = 1
    else:
        zmax = ds[vdim].size
    z_select = pn.widgets.IntSlider(start=0,
                                    end=zmax,
                                    value=1,
                                    name='vertical level',
                                    orientation='vertical',
                                    direction='rtl',
                                    disabled=(vdim is None))

    # bounds for the figures in fraction of the panel,
    # fig_bounds = (0.2, 0.2, 0.8, 0.8)

    @pn.depends(hour_select, z_select)
    def get_quadmesh(hour_select, z_select):
        """
        """
        # if var_varies_vertically:
        #     idx = {'WRFrun': 'control',
        #            'hour': hour_select,
        #            vdim: z_select}
        # else:
        #     idx = {'WRFrun': 'control',
        #            'hour': hour_select}
        if vdim is None:
            qm = ds[varname].sel({'Time': hour_select}).hvplot.quadmesh(
                x='west_east',
                y='south_north',
                z=varname,
                cmap=cmap)
        else:
            qm = ds[varname].sel({'Time': hour_select,
                                  vdim: z_select}).hvplot.quadmesh(
                                      x='west_east',
                                      y='south_north',
                                      z=varname,
                                      cmap=cmap)
        return(qm)

    main_title = '## ' + varname
    the_plot = pn.Column(pn.Row(pn.pane.Markdown(main_title)),
                         pn.Row(get_quadmesh,
                                hour_select,
                                z_select))
    return(the_plot)
