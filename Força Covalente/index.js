import "./d3-force-3d-md.js";

const node_data = [
  { id:0, x: 614.7594242356022, y: 118.22934282922354, r: 4.546777534949488 },
  {id: 1, x: 628.0457626976812, y: 284.00681357087564, r: 4.846936523517396 },
  {id: 2, x: 772.4093106124367, y: 125.39611370240131, r: 6.496580332267502 },
  {id: 3, x: 663.9755973529144, y: 184.28366291900485, r: 4.269674637788155 },
  {id: 4, x: 734.7652403555802, y: 317.6765011625119, r: 17.10323222437136 },
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
//.d3Force("link", d3.forceLink(edges))

// .d3Force('box', () => {
			
//     const CUBE_HALF_SIDE = Graph.nodeRelSize() * N * 0.25;

//     nodes.forEach(node => {
//       const x = node.x || 0, y = node.y || 0, z = node.z || 0;

//       // bounce on box walls
//       if (Math.abs(x) > CUBE_HALF_SIDE) { node.vx *= -1; }
//       if (Math.abs(y) > CUBE_HALF_SIDE) { node.vy *= -1; }
//       if (Math.abs(z) > CUBE_HALF_SIDE) { node.vz *= -1; }
//     });
//   })
  //.nodeAutoColorBy("group")

// Set the data
.graphData({ nodes, links: edges });