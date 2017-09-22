package de.bioforscher.pnator.algorithm;

import de.bioforscher.singa.chemistry.parser.pdb.structures.StructureParser;
import de.bioforscher.singa.chemistry.physical.model.Structure;
import de.bioforscher.singa.core.utility.Resources;
import org.junit.Test;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.UncheckedIOException;
import java.nio.file.Paths;

/**
 * @author cl
 */
public class PNAGeneratorTest {

    private static final Logger logger = LoggerFactory.getLogger(PNAGeneratorTest.class);

    @Test
    public void shouldGeneratePNAForStructures() {

        StructureParser.MultiParser multiParser = StructureParser.online()
                .chainList(Paths.get(Resources.getResourceAsFileLocation("nt_identifier_all.txt")));

        for (int index = 0; index < 20727; index++) {

            try {
                Structure next = multiParser.next();
                logger.info("Converted {}/{} structure {}_{}", index + 1, multiParser.getNumberOfQueuedStructures(),
                        multiParser.getCurrentPdbIdentifier(), multiParser.getCurrentChainIdentifier());
                PNAGenerator.convertToPNAStructure(next);
            } catch (InvalidInputStructure exception) {
                logger.error("Converter failed because an invalid structure. ({}_{})", multiParser.getCurrentPdbIdentifier(), multiParser.getCurrentChainIdentifier());
            } catch (UncheckedIOException exception) {
                logger.error("The PDB identifier {}_{} does not seem to exist", multiParser.getCurrentPdbIdentifier(), multiParser.getCurrentChainIdentifier());
            }

        }
    }

}


