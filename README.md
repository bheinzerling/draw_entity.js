# draw_entity.js

This is a JavaScript package for visualizing entities in text. "Entities" are text spans ("mentions") that are grouped by a certain relation, e.g. the mentions "Paul Allen", "Allen", "Paul" which are all referring to the person Paul Allen in this text:

![Euclidean minimum spanning tree](/images/ex_euclid_mst.png)

There are several options for drawing edges between mentions:

* Euclidean minimum spanning tree: Connects all mentions while minimizing overall line length (pictured above).

* Connect antecedents: connects all mentions in textual order.

* Complete graph: draws all pairwise edges between mentions.

![Connect antecedents](/images/ex_antecedents.png)

![Complete graph](/images/ex_compl_graph.png)

# API:

drawEntity(divIds, style, edgeSelectionFunc)

where

* `divIDs`: The IDs of the divs representing the entity mentions. That is, the divs with these IDs will get connected.

* `style`: The line style used to draw connecting lines, [documented here](https://github.com/sporritt/jsPlumb/wiki/paint-styles)

* `edgeSelectionFunc`: One of `euclideanMST`, `connectAntecedents`, `completeGraph`

## Minimal working example

TODO

## Example

TODO

## Dependencies

* [jsPlumb](https://github.com/sporritt/jsPlumb) (add via `<script>`)

* [jQuery](https://github.com/jquery/jquery) (add via `<script>`)

* [Delaunay](https://github.com/ironwallaby/delaunay) (included)

* [Kruskal](https://github.com/abetusk/kruskal.js) (included)

* [EuclideanMST](https://github.com/abetusk/kruskal.js) (included)
