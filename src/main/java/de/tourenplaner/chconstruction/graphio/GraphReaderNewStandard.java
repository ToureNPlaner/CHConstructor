package de.tourenplaner.chconstruction.graphio;

import de.tourenplaner.chconstruction.graph.RAMGraph;
import fmi.graph.exceptions.NoGraphOpenException;
import fmi.graph.exceptions.NoSuchElementException;
import fmi.graph.maxspeed.Edge;
import fmi.graph.standard.Node;
import fmi.graph.maxspeed.Reader;
import fmi.graph.metaio.MetaData;

import java.io.IOException;
import java.io.InputStream;

/**
 * User: Niklas Schnelle
 * Date: 1/10/13
 * Time: 2:13 PM
 * To change this template use File | Settings | File Templates.
 */
public class GraphReaderNewStandard implements GraphReader {

    /*
    ('motorway','1'),
    ('motorway_link','2'),
    ('trunk','3'),
    ('trunk_link','4'),
    ('primary', '5'),
    ('primary_link','6'),
    ('secondary','7'),
    ('secondary_link','8'),
    ('tertiary','9'),
    ('tertiary_link','10'),
    ('unclassified','11'),
    ('residential','12'),
    ('living_street','13'),
    ('road','14'),
    ('service','15'),
    ('turning_circle','16');
    * */

    private int calcWeight(int length, int roadType, int maxSpeed) {
        double weight;
        if (maxSpeed <= 0){
            switch (roadType) {
                //motorway
                case 1:
                    weight = (length * 1.3) / 1.3;
                    break;
                //motorway link
                case 2:
                    weight = (length * 1.3) / 1.0;
                    break;
                //primary
                case 3:
                    weight = (length * 1.3) / 0.7;
                    break;
                //primary link
                case 4:
                    weight = (length * 1.3) / 0.7;
                    break;
                //secondary
                case 5:
                    weight = (length * 1.3) / 0.65;
                    break;
                //secondary link
                case 6:
                    weight = (length * 1.3) / 0.65;
                    break;
                //tertiary
                case 7:
                    weight = (length * 1.3) / 0.6;
                    break;
                //tertiary link
                case 8:
                    weight = (length * 1.3) / 0.6;
                    break;
                //trunk
                case 9:
                    weight = (length * 1.3) / 0.8;
                    break;
                //trunk link
                case 10:
                    weight = (length * 1.3) / 0.8;
                    break;
                //unclassified
                case 11:
                    weight = (length * 1.3) / 0.25;
                    break;
                //residential
                case 12:
                    weight = (length * 1.3) / 0.45;
                    break;
                //living street
                case 13:
                    weight = (length * 1.3) / 0.3;
                    break;
                //road
                case 14:
                    weight = (length * 1.3) / 0.25;
                    break;
                //service
                case 15:
                    weight = (length * 1.3) / 0.3;
                    break;
                //turning circle
                case 16:
                    weight = (length * 1.3) / 0.3;
                    break;
                default:
                    weight = (length * 1.3) / 0.5;
            }
        } else {
            weight = (length * 1.3) / (((maxSpeed > 130)? 130.0 : (double) maxSpeed)/100.0);
        }


        return (int) weight;
    }

    @Override
    public RAMGraph createRAMGraph(InputStream in) throws IOException {

        Reader r = new Reader();
        MetaData meta = r.read(in);

        Node n;
        Edge e;
        int nodes=0;
        int edges =0;
        int nofNodes;
        int nofEdges;
        long curTime = System.currentTimeMillis();
        try {
            nofNodes = r.getNodeCount();
            nofEdges = r.getEdgeCount();
            RAMGraph graph = new RAMGraph(nofNodes, nofEdges, meta);
            while(r.hasNextNode())
            {
                n = r.nextNode();
                nodes++;
                graph.addNode((float)n.getLat(), (float) n.getLon(), n.getId(), n.getElevation(),(int) n.getOsmId());

                if ((nodes % (nofNodes / 10)) == 0) {
                    System.err.print((10 * nodes / (nofNodes / 10) + "% "));
                }
            }
            System.out.println("Nodes gelesen: "+nodes+" Nodes vorhanden: "+r.getNodeCount());
            while(r.hasNextEdge())
            {
                e = (Edge) r.nextEdge();
                edges++;
                int edgeSource = e.getSource();
                int edgeTarget = e.getTarget();
                int edgeHeight = graph.getHeight(edgeTarget)-graph.getHeight(edgeSource);
                int weight = calcWeight(e.getWeight(), e.getType(), e.getMaxspeed());
                graph.addEdge(e.getSource(), e.getTarget(), weight, e.getWeight(), edgeHeight);

                if ((edges % (nofEdges / 10)) == 0) {
                    System.err.print((10 * edges / (nofEdges / 10) + "% "));
                }
            }
            System.out.println("Edges gelesen: "+edges+" Edges vorhanden: "+r.getEdgeCount());
            System.err.println("Parsing took " + (System.currentTimeMillis() - curTime));
            graph.setupOffsets();

            System.err.println("Read graph with " + nofNodes +
                    " vertices and " + nofEdges + " edges in time " + (System.currentTimeMillis() - curTime) + "ms");
            return graph;
        } catch (NoGraphOpenException ex) {
            throw new IOException(ex);
        } catch (NoSuchElementException ex) {
            throw new IOException(ex);
        }
    }
}
