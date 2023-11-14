import {forceSimulation, forceCollide,forceLink, forceManyBody } from "https://cdn.jsdelivr.net/npm/d3-force@3/+esm";
import { ForceGraph3D } from 'https://cdn.jsdelivr.net/npm/force-graph@2/dist/force-graph.js';

const node_data = [
  { x: 614.7594242356022, y: 118.22934282922354, r: 4.546777534949488 },
  { x: 628.0457626976812, y: 284.00681357087564, r: 4.846936523517396 },
  { x: 772.4093106124367, y: 125.39611370240131, r: 6.496580332267502 },
  { x: 663.9755973529144, y: 184.28366291900485, r: 4.269674637788155 },
  { x: 734.7652403555802, y: 317.6765011625119, r: 17.10323222437136 },
];

const edge_data = [
  { source: 0, target: 1 },
  { source: 1, target: 2 },
  { source: 2, target: 3 },
  { source: 3, target: 4 },
];

const edges = edge_data.map((d) => Object.create(d));

// 1. create a copy of the node data
const nodes = node_data.map((d) => Object.create(d));

// Create the ForceGraph3D instance
const graph = ForceGraph3D()(document.getElementById('graph-container'));

graph
.d3Force("link", forceLink(edges))
.force("charge", forceManyBody().strength(-8))
.force("collide", forceCollide().radius(9))

.d3Force('box', () => {
			
    const CUBE_HALF_SIDE = Graph.nodeRelSize() * N * 0.25;

    nodes.forEach(node => {
      const x = node.x || 0, y = node.y || 0, z = node.z || 0;

      // bounce on box walls
      if (Math.abs(x) > CUBE_HALF_SIDE) { node.vx *= -1; }
      if (Math.abs(y) > CUBE_HALF_SIDE) { node.vy *= -1; }
      if (Math.abs(z) > CUBE_HALF_SIDE) { node.vz *= -1; }
    });
  })
  .nodeAutoColorBy("group")

// Set the data
graphData({ nodes, links: edges });

forceSimulation()
    		.nodes(nodes)