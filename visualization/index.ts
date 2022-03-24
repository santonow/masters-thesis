/**
 * This example shows how to load a GEXF graph file (using the dedicated
 * graphology parser), and display it with some basic map features: Zoom in and
 * out buttons, reset zoom button, and a slider to increase or decrease the
 * quantity of labels displayed on screen.
 */

import Sigma from "sigma";
import { Coordinates, EdgeDisplayData, NodeDisplayData } from "sigma/types";
import Graph from "graphology";

import data from "./data.json";


const graph = new Graph();
graph.import(data)

let sources: Set<string> = new Set<string>();


graph.forEachEdge((edge: string) => {
  const sign: string = graph.getEdgeAttribute(edge, "sign");
  if (sign === "+") {
    graph.setEdgeAttribute(edge, "color", "#99c2ff");
  } else if (sign === "-") {
    graph.setEdgeAttribute(edge, "color", "#ff9999");
  } else {
    graph.setEdgeAttribute(edge, "color", "#bfbfbf");
  }
  graph.getEdgeAttributes(edge).sources.forEach((source: string) => {
    sources.add(source);
  })
})

const buttons = document.getElementById("buttons");
let sourceButtons: Record<string, HTMLButtonElement> = {};
sources.forEach( (source: string) => {
  const button = document.createElement("button") as HTMLButtonElement;
  button.id = source;
  button.textContent = source;
  buttons?.appendChild(button);
  sourceButtons[source] = button;
})


function setClicked(button: HTMLButtonElement) {
  button.style.background = "#d4d4d4";
}

function setUnclicked(button: HTMLButtonElement) {
  button.style.background = "#ffffff";
}


// Retrieve some useful DOM elements:
const container = document.getElementById("sigma-container") as HTMLElement;

const searchInput = document.getElementById("search-input") as HTMLInputElement;
const searchSuggestions = document.getElementById("suggestions") as HTMLDataListElement;

// Instantiate sigma:
const renderer = new Sigma(graph, container, {
  minCameraRatio: 0.1,
  maxCameraRatio: 10,
});

// Type and declare internal state:
interface State {
  hoveredNode?: string;
  searchQuery: string;

  // State derived from query:
  selectedNode?: string;
  suggestions?: Set<string>;

  // State derived from hovered node:
  hoveredNeighbors?: Set<string>;

  // State derived from edge visualization options.
  selectedNetworks?: Set<string>;
}
const state: State = { searchQuery: "", selectedNetworks: new Set<string>(["consensus"])};

// Feed the datalist autocomplete values:
searchSuggestions.innerHTML = graph
  .nodes()
  .map((node) => `<option value="${graph.getNodeAttribute(node, "label")}"></option>`)
  .join("\n");

// Actions:
function setSearchQuery(query: string) {
  state.searchQuery = query;

  if (searchInput.value !== query) searchInput.value = query;

  if (query) {
    const lcQuery = query.toLowerCase();
    const suggestions = graph
      .nodes()
      .map((n) => ({ id: n, label: graph.getNodeAttribute(n, "label") as string }))
      .filter(({ label }) => label.toLowerCase().includes(lcQuery));

    // If we have a single perfect match, them we remove the suggestions, and
    // we consider the user has selected a node through the datalist
    // autocomplete:
    if (suggestions.length === 1 && suggestions[0].label === query) {
      state.selectedNode = suggestions[0].id;
      state.suggestions = undefined;

      // Move the camera to center it on the selected node:
      const nodePosition = renderer.getNodeDisplayData(state.selectedNode) as Coordinates;
      renderer.getCamera().animate(nodePosition, {
        duration: 500,
      });
    }
    // Else, we display the suggestions list:
    else {
      state.selectedNode = undefined;
      state.suggestions = new Set(suggestions.map(({ id }) => id));
    }
  }
  // If the query is empty, then we reset the selectedNode / suggestions state:
  else {
    state.selectedNode = undefined;
    state.suggestions = undefined;
  }

  // Refresh rendering:
  renderer.refresh();
}

function setHoveredNode(node?: string) {
  if (node) {
    state.hoveredNode = node;
    state.hoveredNeighbors = new Set(graph.neighbors(node));
  } else {
    state.hoveredNode = undefined;
    state.hoveredNeighbors = undefined;
  }

  // Refresh rendering:
  renderer.refresh();
}

function setSelectedNetwork(network?: string) {
  let selected_networks: Set<string> = new Set<string>(state.selectedNetworks);
  if (selected_networks.has(network)) {
    selected_networks.delete(network);
  } else {
    selected_networks.add(network);
  }
  state.selectedNetworks = selected_networks;

  renderer.refresh();
}

// Bind search input interactions:
searchInput.addEventListener("input", () => {
  setSearchQuery(searchInput.value || "");
});
searchInput.addEventListener("blur", () => {
  setSearchQuery("");
});

// Bind graph interactions:
renderer.on("enterNode", ({ node }) => {
  setHoveredNode(node);
});
renderer.on("leaveNode", () => {
  setHoveredNode(undefined);
});

function setClickedState(button: HTMLButtonElement, source: string) {
  if (state.selectedNetworks.has(source)) {
    setClicked(button);
  } else {
    setUnclicked(button);
  }
}

Object.entries(sourceButtons).forEach(([source, button]) => {
  setClickedState(button, source);
  button.addEventListener("click", () => {
    setSelectedNetwork(source);
    setClickedState(button, source);
    console.log(`Clicked ${source}`);
    console.log(state.selectedNetworks);
  })
});


function showNeighborNodeLabel(node: string) {
  return state.hoveredNeighbors
      && state.selectedNode
      && state.hoveredNeighbors.has(node)
      && !graph.getEdgeAttribute(graph.edge(state.selectedNode, node), "hidden");
}

// Render nodes accordingly to the internal state:
// 1. If a node is selected, it is highlighted
// 2. If there is query, all non-matching nodes are greyed
// 3. If there is a hovered node, all non-neighbor nodes are greyed
renderer.setSetting("nodeReducer", (node, data) => {
  const res: Partial<NodeDisplayData> = { ...data };

  if (state.hoveredNeighbors && !state.hoveredNeighbors.has(node) && state.hoveredNode !== node) {
    res.label = "";
    res.color = "#f6f6f6";
  }

  if (state.selectedNode === node || showNeighborNodeLabel(node)) {
    res.highlighted = true;
  } else if (state.suggestions && !state.suggestions.has(node)) {
    res.label = "";
    res.color = "#f6f6f6";
  }

  return res;
});

function showEdge(set: Set<string>, subset: Set<string>) {
  if ((subset.has("consensus") && set.has("consensus"))) {
    return true;
  }
  for (let elem of subset.values()) {
    if (!set.has(elem)) {
      return false;
    }
  }
  return true;
}

// Render edges accordingly to the internal state:
// 1. If a node is hovered, the edge is hidden if it is not connected to the
//    node
// 2. If there is a query, the edge is only visible if it connects two
//    suggestions
renderer.setSetting("edgeReducer", (edge, data) => {
  const res: Partial<EdgeDisplayData> = { ...data };

  if (!showEdge(state.selectedNetworks, new Set<string>(graph.getEdgeAttribute(edge, "sources")))) {
    res.hidden = true
  } else if (
      state.hoveredNode && !graph.hasExtremity(edge, state.hoveredNode)
  ) {
    res.hidden = true;
  } else if (
      state.suggestions && (!state.suggestions.has(graph.source(edge)) || !state.suggestions.has(graph.target(edge)))
  ) {
    res.hidden = true;
  }
  return res;
});

const camera = renderer.getCamera();


camera.setState({
    angle: 0.3
})
