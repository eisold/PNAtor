package de.bioforscher.pnator.algorithm;

import de.bioforscher.singa.core.utility.Resources;
import de.bioforscher.singa.structure.model.interfaces.Structure;
import de.bioforscher.singa.structure.parser.pdb.structures.StructureParser;
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

        StructureParser.MultiParser multiParser = StructureParser.pdb()
                .chainList(Paths.get(Resources.getResourceAsFileLocation("nt_identifier_short.txt")));

        for (int index = 0; index < 25187; index++) {

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


