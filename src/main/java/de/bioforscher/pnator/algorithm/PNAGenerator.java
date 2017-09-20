package de.bioforscher.pnator.algorithm;

import de.bioforscher.singa.chemistry.descriptive.elements.Element;
import de.bioforscher.singa.chemistry.parser.pdb.structures.StructureParser;
import de.bioforscher.singa.chemistry.parser.pdb.structures.StructureParserOptions;
import de.bioforscher.singa.chemistry.parser.pdb.structures.StructureWriter;
import de.bioforscher.singa.chemistry.physical.atoms.Atom;
import de.bioforscher.singa.chemistry.physical.atoms.AtomName;
import de.bioforscher.singa.chemistry.physical.atoms.RegularAtom;
import de.bioforscher.singa.chemistry.physical.leaves.Nucleotide;
import de.bioforscher.singa.chemistry.physical.model.Structure;
import de.bioforscher.singa.mathematics.vectors.Vector3D;
import de.bioforscher.singa.mathematics.vectors.Vectors3D;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.IOException;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.List;
import java.util.Optional;

import static de.bioforscher.pnator.algorithm.AtomNameEquivalents.*;
import static de.bioforscher.singa.chemistry.descriptive.elements.ElementProvider.*;

/**
 * @author Alexander Eisold
 * @version 1.0
 */

public class PNAGenerator {

    private static final Logger logger = LoggerFactory.getLogger(PNAGenerator.class);

    private static int addedAtomIndex = 1000;
    private static final double C_O_DOUBLE_BOND_DISTANCE = 1.21;

    private static int backboneFailCount = 0;

    public static void main(String[] args) {

        /*
         * TODO: decide whether structure is dna or rna or hybrid
         * TODO: remove hydrogen atoms if present
         */

        Structure structure = StructureParser.online()
                .pdbIdentifier("1A1V")
                .everything()
                .setOptions(StructureParserOptions.withSettings(StructureParserOptions.Setting.OMIT_HYDROGENS))
                .parse();
        /*
         Structure structure = StructureParser.local().inputStream(Thread.currentThread().getContextClassLoader()
         .getResourceAsStream("structure_examples/example1.pdb")).allModels().parse();
         logger.info("Parsing structure {}.", structure.getPdbIdentifier());
         */
        convertToPNAStructure(structure);

        try {
            StructureWriter.writeBranchSubstructure(structure.getFirstModel(), Paths.get("/tmp/1A1V_PNA.pdb"));
        } catch (IOException e) {
            e.printStackTrace();
        }

        /*
         StructureViewer.colorScheme = ColorScheme.BY_ELEMENT;
         StructureViewer.structure = structure;
         Application.launch(StructureViewer.class);
         */

    }

    public static Structure convertToPNAStructure(Structure structure) {

        structure.getAllChains().forEach(chain -> {

            List<Nucleotide> nucleotides = chain.getNucleotides();
            logger.info("Collected {} nucleotides for chain {}.", nucleotides.size(), chain.getIdentifier());

            if (nucleotides.isEmpty()) {
                logger.info("Chain {} is no nucleosid, skipping.", chain.getIdentifier());
            } else {
                nucleotides.forEach((Nucleotide nucleotide) -> {

                    Optional<Atom> firstPhosphateOptional = FIRST_BACKBONE_PHOSPHATE.getAtomFrom(nucleotide);
                    Optional<Atom> secondPhosphateOptional = SECOND_BACKBONE_PHOSPHATE.getAtomFrom(nucleotide);
                    Optional<Atom> backbonePosphateOptional = BACKBONE_PHOSPHATE.getAtomFrom(nucleotide);

                    Optional<Atom> backboneCarbonOnePrimeOptional = BACKBONE_CARBON_ONE_PRIME.getAtomFrom(nucleotide);
                    Optional<Atom> backboneCarbonTwoPrimeOptional = BACKBONE_CARBON_TWO_PRIME.getAtomFrom(nucleotide);
                    Optional<Atom> backboneCarbonThreePrimeOptional = BACKBONE_CARBON_THREE_PRIME.getAtomFrom(nucleotide);

                    Optional<Atom> backboneOxygenFourPrimeOptional = BACKBONE_OXYGEN_FOUR_PRIME.getAtomFrom(nucleotide);



                    Optional<Atom> backboneOxygenTwoPrimeOptional = BACKBONE_OXYGEN_TWO_PRIME.getAtomFrom(nucleotide);

                    if (backboneCarbonOnePrimeOptional.isPresent() && backboneCarbonTwoPrimeOptional.isPresent() &&
                            backboneCarbonThreePrimeOptional.isPresent()) {

                        Atom oxygenTow = PNAGenerator
                                .calculateMissingAtoms(nucleotide
                                        .getAtomByName(AtomName
                                                .getAtomNameFromString("C1'")), nucleotide
                                        .getAtomByName(AtomName
                                                .getAtomNameFromString("C3'")), nucleotide
                                        .getAtomByName(AtomName
                                                .getAtomNameFromString("C2'")), false, "O7'", OXYGEN);

                        nucleotide.addNode(oxygenTow);
                        nucleotide.addEdgeBetween(nucleotide.getAtomByName(AtomName.getAtomNameFromString("C3'")), oxygenTow);
                    } else {
                        logger.warn("Could not calculate new backbone atome {}.", "C3'");
                    }


                    if (firstPhosphateOptional.isPresent() && secondPhosphateOptional.isPresent() &&
                            backbonePosphateOptional.isPresent()) {

                        Atom oxygenOne = PNAGenerator.calculateMissingAtoms(firstPhosphateOptional.get(),
                                secondPhosphateOptional.get(), backbonePosphateOptional.get(), true,
                                "O1'", OXYGEN);

                        nucleotide.addNode(oxygenOne);
                        nucleotide.addEdgeBetween(backbonePosphateOptional.get(), oxygenOne);

                    } else {
                        backboneFailCount++;
                        if (backboneFailCount == 1) {
                            logger.warn("Could not calculate backbone for nucleotide {}.", nucleotide);
                        } else {
                            logger.warn("Could not calculate backbone more than once.");
                        }
                    }

                    // remove obsolete atoms
                    firstPhosphateOptional.ifPresent(nucleotide::removeNode);
                    secondPhosphateOptional.ifPresent(nucleotide::removeNode);
                    backboneOxygenFourPrimeOptional.ifPresent(nucleotide::removeNode);
                    //chain.removeNode(nucleotide.getAtomByName(AtomName.O4Pr));

                    // remove obsolete RNA specific atoms
                    backboneOxygenTwoPrimeOptional.ifPresent(nucleotide::removeNode);


                    NucleotideValidator validator = new NucleotideValidator(nucleotide.getIdentifier());
                    nucleotide.getAllAtoms().forEach(atom -> convertAtom(atom, validator));
                    if (!validator.isValid()) {
                        logger.warn("The following names were not replaced: {}", validator.getInvalidNames());
                    } else {
                        logger.debug("The nucleotide {} was replaced successfully.", nucleotide.getIdentifier());
                    }

                });
            }
            backboneFailCount = 0;
        });
        return structure;

    }

    private static void convertAtom(Atom an, NucleotideValidator validator) {

        switch (an.getAtomNameString()) {
            case "O5'": {
                String replacement = "N1'";
                replace(an, NITROGEN, replacement);
                validator.validate(replacement);
                break;
            }
            case "C5'": {
                String replacement = "C2'";
                replace(an, CARBON, replacement);
                validator.validate(replacement);
                break;
            }
            case "C4'": {
                String replacement = "C3'";
                replace(an, CARBON, replacement);
                validator.validate(replacement);
                break;
            }
            case "C3'": {
                String replacement = "N4'";
                replace(an, NITROGEN, replacement);
                validator.validate(replacement);
                break;
            }
            case "C2'": {
                String replacement = "C7'";
                replace(an, CARBON, replacement);
                validator.validate(replacement);
                break;
            }
            case "C1'": {
                String replacement = "C8'";
                replace(an, CARBON, replacement);
                validator.validate(replacement);
                break;
            }
            case "O3'": {
                String replacement = "C5'";
                replace(an, CARBON, replacement);
                validator.validate(replacement);
                break;
            }
            case "P": {
                String replacement = "C'";
                replace(an, CARBON, replacement);
                validator.validate(replacement);
                break;
            }
            default:
                break;
        }
    }

    private static void replace(Atom atom, Element replacementElement, String replacementString) {
        logger.trace("Replacing Atom {} through {}.", atom.getAtomNameString(), replacementString);
        atom.setElement(replacementElement);
        atom.setAtomNameString(replacementString);
    }

    /**
     * @param a first Atom for the calculation of the centroid between Atom a and Atom b
     * @param b second Atom for the calculation of the centroid between Atom a and Atom b
     * @param c an Atom that is bounded with the calculated missing Atom
     * @param forward an boolean that defined the direction of the new normalized vector (true is positive and false
     *                negate the vector)
     * @param name an String of the missing atom
     * @param element an Element equals the missing atom
     */
    public static Atom calculateMissingAtoms(Atom a, Atom b, Atom c, Boolean forward, String name, Element element) {

        Vector3D positionA = a.getPosition();
        Vector3D positionB = b.getPosition();

        Vector3D centroid = Vectors3D.getCentroid(Arrays.asList(positionA, positionB));
        logger.trace("Calculated centroid between {} and {} as {}.", a, b, centroid);

        Vector3D missingAtomPosition = null;
        if (forward) {
            missingAtomPosition = c.getPosition().add(centroid.subtract(c.getPosition()).normalize().multiply(C_O_DOUBLE_BOND_DISTANCE));

        } else {
            missingAtomPosition = c.getPosition().add(centroid.subtract(c.getPosition()).normalize().additivelyInvert().multiply(C_O_DOUBLE_BOND_DISTANCE));

        }

        Atom missingAtom = new RegularAtom(addedAtomIndex++, element, name, missingAtomPosition);

        return missingAtom;
    }
}
