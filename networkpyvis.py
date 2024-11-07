import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
from pyvis.network import Network

# Load the CSV file
csv_file = 'output_file.csv'  # Replace with the actual path
df = pd.read_csv(csv_file)

'''# ---- NetworkX Visualization ----
# Create a NetworkX graph
G = nx.Graph()

# Add edges to the NetworkX graph
for _, row in df.iterrows():
    tf = row['TF']  # Assuming 'TF' is the transcription factor column
    target_gene = row['target']  # Assuming 'Target' is the target gene column
    G.add_edge(tf, target_gene)

# Plot the NetworkX graph
plt.figure(figsize=(10, 10))
pos = nx.spring_layout(G, k=0.15)  # Layout for better spacing
nx.draw(G, pos, with_labels=True, node_size=500, node_color='skyblue', font_size=10, font_weight='bold', edge_color='gray')

# Show the plot
plt.title("Transcription Factor - Gene Interaction Network (NetworkX)")
plt.show()'''

# ---- pyvis Interactive Visualization ----
# Create a pyvis Network object
nt = Network(height='600px', width='100%', bgcolor='#ffffff', font_color='black')

# Add nodes and edges to the pyvis network
for _, row in df.iterrows():
    tf = row['TF']
    target_gene = row['target']
    nt.add_node(tf, label=tf)
    nt.add_node(target_gene, label=target_gene)
    nt.add_edge(tf, target_gene)


# Generate and save/show the interactive pyvis network
nt.show('network.html', notebook=False)  # Saves the visualization as an HTML file
