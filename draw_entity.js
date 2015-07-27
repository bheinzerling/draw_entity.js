;

/* 
 * euclideanmst.js by abetusk
 *
 * https://github.com/abetusk/euclideanmst.js
 *
 * GPLv3
 */

var EuclideanMST;
var Delaunay;

  "use strict";

  var EPSILON = 1.0 / 1048576.0;

  function supertriangle(vertices) {
    var xmin = Number.POSITIVE_INFINITY,
        ymin = Number.POSITIVE_INFINITY,
        xmax = Number.NEGATIVE_INFINITY,
        ymax = Number.NEGATIVE_INFINITY,
        i, dx, dy, dmax, xmid, ymid;

    for(i = vertices.length; i--; ) {
      if(vertices[i][0] < xmin) xmin = vertices[i][0];
      if(vertices[i][0] > xmax) xmax = vertices[i][0];
      if(vertices[i][1] < ymin) ymin = vertices[i][1];
      if(vertices[i][1] > ymax) ymax = vertices[i][1];
    }

    dx = xmax - xmin;
    dy = ymax - ymin;
    dmax = Math.max(dx, dy);
    xmid = xmin + dx * 0.5;
    ymid = ymin + dy * 0.5;

    return [
      [xmid - 20 * dmax, ymid -      dmax],
      [xmid            , ymid + 20 * dmax],
      [xmid + 20 * dmax, ymid -      dmax]
    ];
  }

  function circumcircle(vertices, i, j, k) {
    var x1 = vertices[i][0],
        y1 = vertices[i][1],
        x2 = vertices[j][0],
        y2 = vertices[j][1],
        x3 = vertices[k][0],
        y3 = vertices[k][1],
        fabsy1y2 = Math.abs(y1 - y2),
        fabsy2y3 = Math.abs(y2 - y3),
        xc, yc, m1, m2, mx1, mx2, my1, my2, dx, dy;

    /* Check for coincident points */
    if(fabsy1y2 < EPSILON && fabsy2y3 < EPSILON)
      throw new Error("Eek! Coincident points!");

    if(fabsy1y2 < EPSILON) {
      m2  = -((x3 - x2) / (y3 - y2));
      mx2 = (x2 + x3) / 2.0;
      my2 = (y2 + y3) / 2.0;
      xc  = (x2 + x1) / 2.0;
      yc  = m2 * (xc - mx2) + my2;
    }

    else if(fabsy2y3 < EPSILON) {
      m1  = -((x2 - x1) / (y2 - y1));
      mx1 = (x1 + x2) / 2.0;
      my1 = (y1 + y2) / 2.0;
      xc  = (x3 + x2) / 2.0;
      yc  = m1 * (xc - mx1) + my1;
    }

    else {
      m1  = -((x2 - x1) / (y2 - y1));
      m2  = -((x3 - x2) / (y3 - y2));
      mx1 = (x1 + x2) / 2.0;
      mx2 = (x2 + x3) / 2.0;
      my1 = (y1 + y2) / 2.0;
      my2 = (y2 + y3) / 2.0;
      xc  = (m1 * mx1 - m2 * mx2 + my2 - my1) / (m1 - m2);
      yc  = (fabsy1y2 > fabsy2y3) ?
        m1 * (xc - mx1) + my1 :
        m2 * (xc - mx2) + my2;
    }

    dx = x2 - xc;
    dy = y2 - yc;
    return {i: i, j: j, k: k, x: xc, y: yc, r: dx * dx + dy * dy};
  }

  function dedup(edges) {
    var i, j, a, b, m, n;

    for(j = edges.length; j; ) {
      b = edges[--j];
      a = edges[--j];

      for(i = j; i; ) {
        n = edges[--i];
        m = edges[--i];

        if((a === m && b === n) || (a === n && b === m)) {
          edges.splice(j, 2);
          edges.splice(i, 2);
          break;
        }
      }
    }
  }

/*
 * Implementation of Delaunay Triangulation in JavaScript, by ironwallaby
 *
 * https://github.com/ironwallaby/delaunay
 *
 * License: public domain
 */
  Delaunay = {
    triangulate: function(vertices, key) {
      var n = vertices.length,
          i, j, indices, st, open, closed, edges, dx, dy, a, b, c;

      /* Bail if there aren't enough vertices to form any triangles. */
      if(n < 3)
        return [];

      /* Slice out the actual vertices from the passed objects. (Duplicate the
       * array even if we don't, though, since we need to make a supertriangle
       * later on!) */
      vertices = vertices.slice(0);

      if(key)
        for(i = n; i--; )
          vertices[i] = vertices[i][key];

      /* Make an array of indices into the vertex array, sorted by the
       * vertices' x-position. */
      indices = new Array(n);

      for(i = n; i--; )
        indices[i] = i;

      indices.sort(function(i, j) {
        return vertices[j][0] - vertices[i][0];
      });

      /* Next, find the vertices of the supertriangle (which contains all other
       * triangles), and append them onto the end of a (copy of) the vertex
       * array. */
      st = supertriangle(vertices);
      vertices.push(st[0], st[1], st[2]);
      
      /* Initialize the open list (containing the supertriangle and nothing
       * else) and the closed list (which is empty since we havn't processed
       * any triangles yet). */
      open   = [circumcircle(vertices, n + 0, n + 1, n + 2)];
      closed = [];
      edges  = [];

      /* Incrementally add each vertex to the mesh. */
      for(i = indices.length; i--; edges.length = 0) {
        c = indices[i];

        /* For each open triangle, check to see if the current point is
         * inside it's circumcircle. If it is, remove the triangle and add
         * it's edges to an edge list. */
        for(j = open.length; j--; ) {
          /* If this point is to the right of this triangle's circumcircle,
           * then this triangle should never get checked again. Remove it
           * from the open list, add it to the closed list, and skip. */
          dx = vertices[c][0] - open[j].x;
          if(dx > 0.0 && dx * dx > open[j].r) {
            closed.push(open[j]);
            open.splice(j, 1);
            continue;
          }

          /* If we're outside the circumcircle, skip this triangle. */
          dy = vertices[c][1] - open[j].y;
          if(dx * dx + dy * dy - open[j].r > EPSILON)
            continue;

          /* Remove the triangle and add it's edges to the edge list. */
          edges.push(
            open[j].i, open[j].j,
            open[j].j, open[j].k,
            open[j].k, open[j].i
          );
          open.splice(j, 1);
        }

        /* Remove any doubled edges. */
        dedup(edges);

        /* Add a new triangle for each edge. */
        for(j = edges.length; j; ) {
          b = edges[--j];
          a = edges[--j];
          open.push(circumcircle(vertices, a, b, c));
        }
      }

      /* Copy any remaining open triangles to the closed list, and then
       * remove any triangles that share a vertex with the supertriangle,
       * building a list of triplets that represent triangles. */
      for(i = open.length; i--; )
        closed.push(open[i]);
      open.length = 0;

      for(i = closed.length; i--; )
        if(closed[i].i < n && closed[i].j < n && closed[i].k < n)
          open.push(closed[i].i, closed[i].j, closed[i].k);

      /* Yay, we're done! */
      return open;
    },
    contains: function(tri, p) {
      /* Bounding box test first, for quick rejections. */
      if((p[0] < tri[0][0] && p[0] < tri[1][0] && p[0] < tri[2][0]) ||
         (p[0] > tri[0][0] && p[0] > tri[1][0] && p[0] > tri[2][0]) ||
         (p[1] < tri[0][1] && p[1] < tri[1][1] && p[1] < tri[2][1]) ||
         (p[1] > tri[0][1] && p[1] > tri[1][1] && p[1] > tri[2][1]))
        return null;

      var a = tri[1][0] - tri[0][0],
          b = tri[2][0] - tri[0][0],
          c = tri[1][1] - tri[0][1],
          d = tri[2][1] - tri[0][1],
          i = a * d - b * c;

      /* Degenerate tri. */
      if(i === 0.0)
        return null;

      var u = (d * (p[0] - tri[0][0]) - b * (p[1] - tri[0][1])) / i,
          v = (a * (p[1] - tri[0][1]) - c * (p[0] - tri[0][0])) / i;

      /* If we're outside the tri, fail. */
      if(u < 0.0 || v < 0.0 || (u + v) > 1.0)
        return null;

      return [u, v];
    }
  };


function UnionFind(count) {
  this.roots = new Array(count);
  this.ranks = new Array(count);
  
  for(var i=0; i<count; ++i) {
    this.roots[i] = i;
    this.ranks[i] = 0;
  }
}

var proto = UnionFind.prototype

Object.defineProperty(proto, "length", {
  "get": function() {
    return this.roots.length
  }
})

proto.makeSet = function() {
  var n = this.roots.length;
  this.roots.push(n);
  this.ranks.push(0);
  return n;
}

proto.find = function(x) {
  var roots = this.roots;
  while(roots[x] !== x) {
    var y = roots[x];
    roots[x] = roots[y];
    x = y;
  }
  return x;
}

proto.link = function(x, y) {
  var xr = this.find(x)
    , yr = this.find(y);
  if(xr === yr) {
    return;
  }
  var ranks = this.ranks
    , roots = this.roots
    , xd    = ranks[xr]
    , yd    = ranks[yr];
  if(xd < yd) {
    roots[xr] = yr;
  } else if(yd < xd) {
    roots[yr] = xr;
  } else {
    roots[yr] = xr;
    ++ranks[xr];
  }
}

/* 
 * Implementation of Kruskal's algorithm by abetusk.
 *
 * https://github.com/abetusk/kruskal.js
 *
 * License: GPLv3
 */
var Kruskal;
  "use strict";

  // vertices hold data that will be used in the distance 'metric' function
  // edges holds position in vertices list
  //
  Kruskal = {
    kruskal: function( vertices, edges, metric )
    {
      var set = {};

      var finalEdge = [];

      var forest = new UnionFind( vertices.length );

      var edgeDist = [];
      for (var ind in edges)
      {
        var u = edges[ind][0];
        var v = edges[ind][1];
        var e = { edge: edges[ind], weight: metric( vertices[u], vertices[v] ) };
        edgeDist.push(e);
      }

      edgeDist.sort( function(a, b) { return a.weight- b.weight; } );

      for (var i=0; i<edgeDist.length; i++)
      {
        var u = edgeDist[i].edge[0];
        var v = edgeDist[i].edge[1];

        if ( forest.find(u) != forest.find(v) )
        {
          finalEdge.push( [ u, v ] );
          forest.link( u, v );
        }
      }

      return finalEdge;

    }
  }

  EuclideanMST = {
    euclideanMST : function ( vertices, metric )
    {

      var tri = Delaunay.triangulate(vertices);

      var edges = [];
      var edge_seen = {};

      for (var i=0; i<tri.length; i+=3)
      {
        var e = [ tri[i], tri[i+1], tri[i+2] ];
        e.sort();

        var a = e[0];
        var b = e[1];
        var c = e[2];

        var key = "" + a + "," + b;
        if ( !(key in edge_seen) )
          edges.push( [a,b] );
        edge_seen[key] = 1;

        key = "" + b + "," + c;
        if ( !(key in edge_seen) )
          edges.push( [b,c] );
        edge_seen[key] = 1;

        key = "" + a + "," + c;
        if ( !(key in edge_seen) )
          edges.push( [a,c] );
        edge_seen[key] = 1;

      }

      var mst = 
        Kruskal.kruskal( 
            vertices, 
            edges, 
            metric
        );

      return mst;

    }

  };

// end imports

/*
 * Draw edges connecting the supplied divIds using the supplied line style.
 *
 * The line style is an object as described in the jsPlumb documentation:
 * http://www.jsplumb.org/doc/paint-styles.html
 *
 * E.g., to draw green, dashed lines 3 pixel wide, use something like this:
 * {strokeStyle: '#2ECC40', dashstyle: "2 3", lineWidth:3}
 *
 * How the edges are selected is specified using edgeSelectionFunc, which is
 * one of: 
 *
 *      euclideanMST: edges that connect the divIds form a Euclidiean minimum spannning tree
 *      completeGraph: draw all possible edges between the divIds
 *      connectAntecedents: draw edges connect each div ID to its antecedent div ID
 */
function drawEntity(divIds, style, edgeSelectionFunc) {
    var edges = edgeSelectionFunc(divIds);
    drawEdges(edges, style);
}

/*
 * Given a list of entities in the form of a list of lists of divIDs, draw each entity
 * using the supplied line style and edgeSelectionFunc, as described in the drawEntity
 * documentation.
 */
function drawEntities(divIdsList, style, edgeSelectionFunc) {
    divIdsList.forEach(function (divIds) {
        drawEntity(divIds, style, edgeSelectionFunc);
    });
}

/*
 * Computes the Euclidan minimum spanning tree of the vertices represented by the supplied div IDs.
 *
 * The MST is returned as a list of edges, where an edge is a pair of div IDs,
 */
function euclideanMST(divIds) {
    var vertices = positionsOfDivs(divIds);
    var edgeIdxPairs = EuclideanMST.euclideanMST(vertices, squaredEuclideanDistance);        
    // eucldiean;ST seems to erroneously return an empty tree for complete 2-node graphs
    if (edgeIdxPairs.length == 0 && divIds.length == 2) {
        edgeIdxPairs = [[0, 1]];
    }
    return edgeIdxPairs.map(function(edgeIdxPair) {return [divIds[edgeIdxPair[0]], divIds[edgeIdxPair[1]]];});
}

/*
 * Computes the complete graph of the vertices represented by the supplied div IDs.
 *
 * The complete graph is returned as a list of edges, where an edge is a pair of div IDs,
 */
function completeGraph(divIds) {
    // http://codereview.stackexchange.com/questions/75658/pairwise-combinations-of-an-array-in-javascript
    var pairs = new Array((divIds.length * (divIds.length - 1)) / 2),
       pos = 0;

   for (var i = 0; i < divIds.length; i++) {
       for (var j = i + 1; j < divIds.length; j++) {
           pairs[pos++] = [divIds[i], divIds[j]];
       }
   }
   return pairs;
}

/*
 * Computes a graph that connects the vertices represented by the supplied div IDs according to their textual order.
 *
 * The graph is returned as a list of edges, where an edge is a pair of div IDs,
 */
function connectAntecedents(divIds) {
    divIds.sort(function(divId1, divId2) {
        var pos1 = positionOfDiv(divId1);
        var pos2 = positionOfDiv(divId2);
        if (pos1[0] == pos2[0]) {
            return pos1[1] - pos2[1];
        }
        return pos1[0] - pos2[0];
    });
    return zip(divIds.slice(0, -1), divIds.slice(1));
}

function positionOfDiv(divId) {
    var offset = $("#" + divId).offset();
    return [offset.top, offset.left];
}

function positionsOfDivs(divIds) {
    return divIds.map(positionOfDiv);
}

function squaredEuclideanDistance(a,b) {
  return (a[0]-b[0])*(a[0]-b[0]) + (a[1]-b[1])*(a[1]-b[1]); 
}

function squaredDivDistance(divId1, divId2) {
    return squaredEuclideanDistance(positionOfDiv(divId1), positionOfDiv(divId2));
}

function zip(a, b) {
    return a.map(function (e, i) {
        return [a[i], b[i]];
    });
}

function anchorsForEdge(edge) {
    positionSource = positionOfDiv(edge[0]);
    positionTarget = positionOfDiv(edge[1]);
    if (positionSource[0] > positionTarget[0] ) {
        return ['Bottom', 'Top'];
    }
    if (positionSource[0] < positionTarget[0]) {
        return ['Top', 'Bottom'];
    }
    return ['Top', 'Top'];
}

function drawEdge(edge, style) {
    var anchors = anchorsForEdge(edge);
    jsPlumb.connect({
        source: edge[0],
        target: edge[1],
        anchor: anchors, paintStyle: style, endpoint: "Blank", connector: ['Bezier', {curviness:10}],
        overlays:[["Label", {cssClass: "label", label: ""}]]
        });
}

function drawEdges(edges, style) {
    edges.forEach(function (edge) {
        drawEdge(edge, style);
    });
}

function findClosestPair(divIds1, divIds2) {
    var minimum = 10e9;
    var currentPair = [];
    divIds1.forEach(function (divId1) { 
        divIds2.forEach(function (divId2) {
            distance = squaredDivDistance(divId1, divId2);
            if (distance < minimum) {
                minimum = distance;
                currentPair = [divId1, divId2];
            }
        });
    });
    return currentPair;
}

