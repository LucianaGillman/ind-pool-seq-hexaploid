import plotly.express as px
import pandas as pd

mds_data = pd.read_csv('MDStable_pool.csv')

fig = px.scatter_3d(mds_data, x='Factor1', y='Factor2', z='Factor3',
                    color='Accession',symbol="Accession",
                    labels={"Factor1": "Factor 1 (22%)",
                     "Factor2": "Factor 2 (17%)", "Factor3": "Factor 3 (11%)"},
                    range_x = [0.2,-0.2],range_y = [0.2,-0.2],range_z = [0.2,-0.2],
			width=1200, height=800)      
fig.update_layout(
   		  scene=dict(
             		    xaxis=dict(
				backgroundcolor="lightgoldenrodyellow",
 				gridcolor="blue",
                                showbackground=True,
                                zerolinecolor="black",),
                            yaxis = dict(
                                backgroundcolor="white",
                                gridcolor="blue",
                                showbackground=True,
                                zerolinecolor="black"),
                            zaxis = dict(
                                backgroundcolor="lightcyan",
                                gridcolor="blue",
                                showbackground=True,
                                zerolinecolor="black",),),
                     )         				 

import kaleido

fig.write_image("MDS_ind_axes_pool.png", engine= "kaleido")        				 

fig.write_image("MDS_ind_axes_pool.pdf", engine= "kaleido")
fig.write_image("MDS_ind_axes_pool.jpg", engine= "kaleido")
fig.write_image("MDS_ind_axes_pool.eps", engine= "kaleido")
