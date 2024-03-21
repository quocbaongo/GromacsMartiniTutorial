from jupyter_dash import JupyterDash
import dash
from dash.dependencies import Input, Output, State
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
#import ipywidgets as widgets
from pathlib import Path
import json
import subprocess
import argparse
import os 
from scipy.spatial import distance

def validate_file(f):
	if not os.path.exists(f):
		# Argparse uses the ArgumentTypeError to give a rejection message like:
		# error: argument input: x does not exist
		raise argparse.ArgumentTypeError("{0} does not exist".format(f))
	return f


if __name__ == "__main__":

	# flag list
	parser = argparse.ArgumentParser()
	parser = argparse.ArgumentParser(description = 'Interacting with the principle component map to extract favorable states')
	parser.add_argument("--input", type=validate_file, help='Input file (2dproj.xvg)')	
	parser.add_argument('--traj', help='Trajectory (traj.xtc)')
	parser.add_argument('--tpr', help='Topology (topol.tpr)')
	parser.add_argument('--dt', type=int, help='Time step (ps)')
	
	args = parser.parse_args()


	# Raw 2d PCA file
	PCAFile=args.input

	# Open PCA file
	with open(PCAFile, 'r') as f:  
		PCAFileContent = f.readlines() 	

	
	# Reading PCAFile content
	PCAData = []
	
	for line in PCAFileContent:
		if not (line.startswith('#') or line.startswith('@')):
			temp = [float(line.split()[0]), float(line.split()[1])]
			PCAData.append(temp)
		
	print(len(PCAData))

	# Shorten the list	
	x = [i[0] for i in PCAData]
	y = [i[1] for i in PCAData]
	# Creating bins
	x_min = np.min([i[0] for i in PCAData])
	x_max = np.max([i[0] for i in PCAData])
  
	y_min = np.min([i[1] for i in PCAData])
	y_max = np.max([i[1] for i in PCAData])
  
	x_bins = np.linspace(x_min, x_max, 100)
	y_bins = np.linspace(y_min, y_max, 100)
	

	H, xedges, yedges = np.histogram2d(x, y, bins=(x_bins, y_bins))
	H = H.T

	# Plotting
	# Plotting energy landscape
	fig = go.Figure(data=[go.Surface(z=H, x=x_bins, y=y_bins)])
	fig.update_layout(title='Energy landscape', autosize=False,
			width=1000, height=1000,
			margin=dict(l=65, r=50, b=65, t=90),
			scene = dict(xaxis_title='PC1',
			yaxis_title='PC2'))


	# Working with Dash
	clicked = []

	# Build App
	app = JupyterDash(__name__)
	app.layout = dash.html.Div(
				[dash.dcc.Graph(id="fig",
						figure=fig,
						),
				dash.html.Div(id="debug"),
				]
				)


	@app.callback(
		Output("debug", "children"),
		Input("fig", "clickData"),
			)
			
	def point_clicked(clickData):
		global clicked
		clicked.append(clickData)

		return json.dumps(clickData)

	# Run app and display result inline in the notebook
	app.run_server(mode="inline",host='127.0.0.1')


	# Start extracting structure from click
	# PC1 data from 2dproj.xvg file
	PC1 = x
	PC2 = y
	commands = []	
	PDBOut = []

	for i in range(len(clicked)):
		if clicked[i] != None:
		
			#print(clicked[i]['points'][0]['x'],clicked[i]['points'][0]['y'],clicked[i]['points'][0]['z'])
			distances = [distance.euclidean([clicked[i]['points'][0]['x'],clicked[i]['points'][0]['y']], 
					[pc1, pc2]) for pc1, pc2 in zip(PC1, PC2)]		
			
	
			index = distances.index(min(distances))
			indexToTime = index * args.dt
			#print(index)
			#print(indexToTime)	

			#print(indexToTime, PC1[index], PC2[index])			



			# Extracting using gromacs
			p = subprocess.Popen(['gmx_mpi', 'trjconv', 
						'-f', args.traj, '-s', args.tpr,
						'-dump', str(indexToTime), '-o', 'frame_' + str(PC1[index]) + '_' + str(PC2[index]) + '_' + str(indexToTime) + '.pdb'],
					stdin=subprocess.PIPE)		
			p.communicate(b'1\n')
			p.wait()

			commands.append('Command just used: ' + 'gmx_mpi trjconv -f ' + args.traj + ' -s '+ args.tpr + ' -dump ' + str(indexToTime) + ' -o frame_' + str(PC1[index]) + '_' + str(PC2[index]) + '_' + str(indexToTime) + '.pdb')		
			
			# Put extracted pdb file in a list
			PDBOut.append('frame_' + str(PC1[index]) + '_' + str(PC2[index]) + 
					'_' + str(indexToTime) + '.pdb')
		
	print(clicked)
	print()
	for i in commands:
		print(i)	








