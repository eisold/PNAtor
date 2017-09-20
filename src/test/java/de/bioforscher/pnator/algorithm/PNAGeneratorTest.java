package de.bioforscher.pnator.algorithm;

import de.bioforscher.singa.chemistry.parser.pdb.structures.StructureParser;
import de.bioforscher.singa.chemistry.physical.model.Structure;
import de.bioforscher.singa.core.utility.Resources;
import org.junit.Test;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.nio.file.Paths;

/**
 * @author cl
 */
public class PNAGeneratorTest {

    private static final Logger logger = LoggerFactory.getLogger(PNAGeneratorTest.class);

    @Test
    public void shouldGeneratePNAForStructures() {

        StructureParser.MultiParser multiParser = StructureParser.online()
                .chainList(Paths.get(Resources.getResourceAsFileLocation("nt_identifier.txt")));

        for (int index = 0; index < 100; index++) {
            Structure next = multiParser.next();
            try {
                logger.info("Converted {}/{} structure {}", index,multiParser.getNumberOfQueuedStructures(),
                        multiParser.getCurrentPdbIdentifier());
                PNAGenerator.convertToPNAStructure(next);
            } catch (InvalidInputStructure exception) {
                logger.error("Converter failed because an invalid structure. ({})", next.getPdbIdentifier());
            }
        }

    }


}