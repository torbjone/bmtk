# Allen Institute Software License - This software license is the 2-clause BSD license plus clause a third
# clause that prohibits redistribution for commercial purposes without further permission.
#
# Copyright 2017. Allen Institute. All rights reserved.
#
# Redistribution and use in source and binary forms, with or without modification, are permitted provided that the
# following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following
# disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following
# disclaimer in the documentation and/or other materials provided with the distribution.
#
# 3. Redistributions for commercial purposes are not permitted without the Allen Institute's written permission. For
# purposes of this license, commercial purposes is the incorporation of the Allen Institute's software into anything for
# which you will charge fees or other compensation. Contact terms@alleninstitute.org for commercial licensing
# opportunities.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,
# INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
# WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
import os
import numpy as np
import types
import csv

from node_pool import NodePool
from connection_map import ConnectionMap
from node_set import NodeSet
from id_generator import IDGenerator


class Network (object):
    def __init__(self, name, **network_props):
        if len(name) == 0:
            raise Exception('Network name missing.')

        self._network_name = name

        self._nnodes = 0
        self._nodes_built = False
        self._nedges = 0
        self._edges_built = False
        
        self._node_sets = []
        self.__external_node_sets = []
        self.__node_id_counter = 0

        self._node_types_properties = {}
        self._node_types_columns = set(['node_type_id'])
        # self._edge_type_properties = {}
        # self._edge_types_columns = set(['edge_type_id'])
        self._connection_maps = []
        #self._connection_maps = ConnectionTable()

        self._node_id_gen = IDGenerator()
        self._node_type_id_gen = IDGenerator(100)
        self._edge_type_id_gen = IDGenerator(100)

        #self._connection_table = []
        #self._source_networks = []
        #self._target_networks = []
        self._network_conns = set()

    @property
    def name(self):
        return self._network_name

    @property
    def nodes_built(self):
        return self._nodes_built

    @property
    def edges_built(self):
        return self._edges_built

    @property
    def nnodes(self):
        raise NotImplementedError

    @property
    def nedges(self):
        raise NotImplementedError

    def get_connections(self):
        return self._connection_maps

    def _add_node_type(self, props):
        node_type_id = props.get('node_type_id', None)
        if node_type_id is None:
            node_type_id = self._node_type_id_gen.next()
        else:
            if node_type_id in self._node_types_properties:
                raise Exception('node_type_id {} already exists.'.format(node_type_id))
            self._node_type_id_gen.remove_id(node_type_id)

        props['node_type_id'] = node_type_id
        self._node_types_properties[node_type_id] = props

    def add_nodes(self, N=1, **properties):
        self._clear()

        # categorize properties as either a node-params (for nodes file) or node-type-property (for node_types files)
        node_params = {}
        node_properties = {}
        for prop_name, prop_value in properties.items():
            if isinstance(prop_value, (list, np.ndarray)):  # TODO: what about pandas series
                n_props = len(prop_value)
                if n_props != N:
                    raise Exception('Trying to pass in array of length {} into N={} nodes'.format(n_props, N))
                node_params[prop_name] = prop_value

            elif isinstance(prop_value, types.GeneratorType):
                vals = list(prop_value)
                assert(len(vals) == N)
                node_params[prop_name] = vals

            else:
                node_properties[prop_name] = prop_value
                self._node_types_columns.add(prop_name)

        # If node-type-id exists, make sure there is no clash, otherwise generate a new id.
        if 'node_type_id' in node_params:
            raise Exception('There can be only one "node_type_id" per set of nodes.')

        self._add_node_type(node_properties)
        self._node_sets.append(NodeSet(N, node_params, node_properties))

    def add_edges(self, source=None, target=None, connection_rule=1, connection_params=None, iterator='one_to_one',
                  **edge_type_properties):
        # TODO: check edge_type_properties for 'edge_type_id' and make sure there isn't a collision. Otherwise create
        #       a new id.
        if not isinstance(source, NodePool):
            source = NodePool(self, **source or {})

        if not isinstance(target, NodePool):
            target = NodePool(self, **target or {})

        self._network_conns.add((source.network_name, target.network_name))

        # TODO: make sure that they don't add a dictionary or some other wried property type.
        edge_type_id = edge_type_properties.get('edge_type_id', None)
        if edge_type_id is None:
            edge_type_id = self._edge_type_id_gen.next()
            edge_type_properties['edge_type_id'] = edge_type_id
        elif edge_type_id in self._edge_type_id_gen:
            raise Exception('edge_type_id {} already exists.'.format(edge_type_id))
        else:
            self._edge_type_id_gen.remove_id(edge_type_id)

        edge_type_properties['source_query'] = source.filter_str
        edge_type_properties['target_query'] = target.filter_str

        # self._edge_types_columns.update(edge_type_properties.keys())
        connection = ConnectionMap(source, target, connection_rule, connection_params, iterator, edge_type_properties)
        self._connection_maps.append(connection)
        # self._connection_maps.add(source.network_name, target.network_name, connection)
        return connection

    def nodes(self, **properties):
        if not self.nodes_built:
            self._build_nodes()

        return NodePool(self, **properties)

    def edges(self, targets=None, sources=None):
        if not self.edges_built:
            self.build()

        return list(self._edges_iter(targets=targets))

    def clear(self):
        self._nodes_built = False
        self._edges_built = False
        self._clear()

    def _node_id(self, N):
        for i in xrange(N):
            yield self.__node_id_counter
            self.__node_id_counter += 1

    def _build_nodes(self):
        """Builds or rebuilds all the nodes, clear out both node and edge sets."""
        # print 'build_nodes'
        self._clear()
        self._initialize()

        for ns in self._node_sets:
            nodes = ns.build(nid_generator=self._node_id)
            self._add_nodes(nodes)
        self._nodes_built = True

    def __build_edges(self):
        """Builds network edges"""
        if not self.nodes_built:
            # only rebuild nodes if necessary.
            self._build_nodes()

        for i, conn_map in enumerate(self._connection_maps):
            # print conn_map
            self._add_edges(conn_map, i)

        self._edges_built = True

    def build(self, force=False):
        """ Builds nodes (assigns gids) and edges.

        Args:
            force (bool): set true to force complete rebuilding of nodes and edges, if nodes() or save_nodes() has been
                called before then forcing a rebuild may change gids of each node.
        """

        # if nodes() or save_nodes() is called by user prior to calling build() - make sure the nodes
        # are completely rebuilt (unless a node set has been added).
        if force:
            self._clear()
            self._initialize()
            self._build_nodes()

        # always build the edges.
        self.__build_edges()

    def save_nodes(self, nodes_file_name, node_types_file_name):
        raise NotImplementedError()

    def import_nodes(self, nodes_file_name, node_types_file_name):
        raise NotImplementedError()

    def save_edges(self, edges_file_name=None, edge_types_file_name=None, output_dir='.', src_network=None,
                   trg_network=None, force_build=True, force_overwrite=False):
        # Make sure edges exists and are built
        if len(self._connection_maps) == 0:
            print("Warning: no edges have been made for this network, skipping saving.")
            return

        if self._edges_built is False:
            if force_build:
                print("Message: building edges")
                self.__build_edges()
            else:
                print("Warning: Edges are not built. Either call build() or use force_build parameter. Skip saving.")
                return

        network_params = [(s, t, s+'_'+t+'_edges.h5', s+'_'+t+'_edge_types.csv') for s, t in list(self._network_conns)]
        if src_network is not None:
            network_params = [p for p in network_params if p[0] == src_network]

        if trg_network is not None:
            network_params = [p for p in network_params if p[1] == trg_network]

        if len(network_params) == 0:
            print("Warning: couldn't find connections. Skip saving.")
            return

        if (edges_file_name or edge_types_file_name) is not None:
            network_params = [(network_params[0][0], network_params[0][1], edges_file_name, edge_types_file_name)]

        if not os.path.exists(output_dir):
            os.mkdir(output_dir)

        for p in network_params:
            if p[3] is not None:
                self._save_edge_types(os.path.join(output_dir, p[3]), p[0], p[1])

            if p[2] is not None:
                self._save_edges(os.path.join(output_dir, p[2]), p[0], p[1])

    def _save_edge_types(self, edge_types_file_name, src_network, trg_network):

        # Get edge-type properties for connections with matching source/target networks
        matching_et = [c.edge_type_properties for c in self._connection_maps
                       if c.source_network_name == src_network and c.target_network_name == trg_network]

        # Get edge-type properties that are only relevant for this source-target network pair
        cols = ['edge_type_id', 'target_query', 'source_query']  #  manditory and should come first
        merged_keys = [k for et in matching_et for k in et.keys() if k not in cols]
        cols += list(set(merged_keys))

        # Write to csv
        with open(edge_types_file_name, 'w') as csvfile:
            csvw = csv.writer(csvfile, delimiter=' ')
            csvw.writerow(cols)
            for edge_type in matching_et:
                csvw.writerow([edge_type.get(cname, 'NULL') for cname in cols])

    def _save_edges(self, edges_file_name, src_network, trg_network):
        raise NotImplementedError

    def _initialize(self):
        raise NotImplementedError

    def _add_nodes(self, node_tuples):
        raise NotImplementedError

    def _add_edges(self, edge_tuples, i):
        raise NotImplementedError

    def _clear(self):
        raise NotImplementedError

    def _nodes_iter(self, nids=None):
        raise NotImplementedError

    def _edges_iter(targets=None, sources=None):
        raise NotImplementedError


class ConnectionTable(object):
    def __init__(self):
        self.__targets = {}
        self.__sources = {}
        self.__connections = []

    def add(self, source_network, target_network, connection_map):
        # TODO: If the source/target are network objects we can get the network_name
        assert(isinstance(source_network, basestring))
        assert(isinstance(target_network, basestring))
        assert(isinstance(connection_map, ConnectionMap))

        if source_network not in self.__sources:
            self.__sources[source_network] = []
        if target_network not in self.__targets:
            self.__targets[target_network] = []

        cm_index = len(self.__connections)
        self.__connections.append(connection_map)
        self.__sources[source_network].append(cm_index)
        self.__targets[target_network].append(cm_index)

    def get(self, source_network=None, target_network=None):
        # TODO: Add warning if source/target network is not found
        cm_indicies = set(range(len(self.__connections)))
        if source_network is not None:
            cm_indicies &= set(self.__sources.get(source_network, []))

        if target_network is not None:
            cm_indicies &= set(self.__targets.get(target_network, []))

        return self.__connections[cm_indicies]





