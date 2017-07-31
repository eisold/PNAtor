package de.bioforscher.pnator.algorithm;

import de.bioforscher.singa.chemistry.parser.pdb.structures.StructureParser;
import de.bioforscher.singa.chemistry.physical.model.Structure;


public class PNAGenerator {

        Structure structure = StructureParser.online()
                .pdbIdentifier("1BNA")
                .parse();

}