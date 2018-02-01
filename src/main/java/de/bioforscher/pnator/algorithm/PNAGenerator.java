package de.bioforscher.pnator.algorithm;

import de.bioforscher.singa.mathematics.vectors.Vector3D;
import de.bioforscher.singa.mathematics.vectors.Vectors3D;
import de.bioforscher.singa.structure.elements.Element;
import de.bioforscher.singa.structure.model.identifiers.LeafIdentifier;
import de.bioforscher.singa.structure.model.interfaces.Atom;
import de.bioforscher.singa.structure.model.interfaces.LeafSubstructure;
import de.bioforscher.singa.structure.model.interfaces.Nucleotide;
import de.bioforscher.singa.structure.model.interfaces.Structure;
import de.bioforscher.singa.structure.model.oak.OakAtom;
import de.bioforscher.singa.structure.model.oak.OakNucleotide;
import de.bioforscher.singa.structure.model.oak.OakStructure;
import de.bioforscher.singa.structure.parser.pdb.structures.StructureParser;
import de.bioforscher.singa.structure.parser.pdb.structures.StructureWriter;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.IOException;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.List;
import java.util.Optional;

import static de.bioforscher.pnator.algorithm.AtomNameEquivalents.*;
import static de.bioforscher.singa.structure.elements.ElementProvider.CARBON;
import static de.bioforscher.singa.structure.elements.ElementProvider.NITROGEN;
import static de.bioforscher.singa.structure.elements.ElementProvider.OXYGEN;

/**
 * @author Alexander Eisold
 * @version 1.0
 */

public class PNAGenerator {

    private static final Logger logger = LoggerFactory.getLogger(PNAGenerator.class);

    private static int addedAtomIndex = 0;
    private static final double C_O_DOUBLE_BOND_DISTANCE = 1.21;

    private static int backboneFailCount = 0;

    public static void main(String[] args) {

        /*
         * TODO: decide whether structure is dna or rna or hybrid
         * TODO: remove hydrogen atoms if present
         */
/*
        Structure structure = StructureParser.online()
                .pdbIdentifier("1OLD")
                .everything()
                .setOptions(StructureParserOptions.withSettings(StructureParserOptions.Setting.OMIT_HYDROGENS))
                .parse();
        */
        OakStructure structure = (OakStructure) StructureParser.local().inputStream(Thread.currentThread().getContextClassLoader()
                .getResourceAsStream("structure_examples/model_2501_only_aptamer.pdb")).allModels().parse();
        logger.info("Parsing structure {}.", structure.getPdbIdentifier());
        addedAtomIndex = structure.getLastAddedAtomIdentifier();
        convertToPNAStructure(structure);

        try {
            StructureWriter.writeLeafSubstructureContainer(structure.getFirstModel(), Paths.get("/home/aeisold/Documents/work/PNA/PNAtor/output/model_2501_only_aptamer_PNA.pdb"));
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

            List<Nucleotide> nucleotides = chain.getAllNucleotides();
            logger.info("Collected {} nucleotides for chain {}.", nucleotides.size(), chain.getChainIdentifier());

            if (nucleotides.isEmpty()) {
                logger.info("Chain {} is no nucleosid, skipping.", chain.getChainIdentifier());
            } else {
                nucleotides.stream().map(OakNucleotide.class::cast)
                        .forEach(nucleotide -> {

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

                                OakAtom oxygenTwo = PNAGenerator
                                        .calculateMissingAtoms(nucleotide.getAtomByName("C1'").get(),
                                                nucleotide.getAtomByName("C3'").get(),
                                                nucleotide.getAtomByName("C2'").get(), false, "O7'", OXYGEN);

                                nucleotide.addAtom(oxygenTwo);
                                nucleotide.addBondBetween((OakAtom) nucleotide.getAtomByName("C3'").get(), oxygenTwo);
                            } else {
                                logger.warn("Could not calculate new backbone atome {}.", "C3'");
                            }


                            if (firstPhosphateOptional.isPresent() && secondPhosphateOptional.isPresent() &&
                                    backbonePosphateOptional.isPresent()) {

                                OakAtom oxygenOne = PNAGenerator.calculateMissingAtoms(firstPhosphateOptional.get(),
                                        secondPhosphateOptional.get(), backbonePosphateOptional.get(), true,
                                        "O1'", OXYGEN);


                                LeafIdentifier leafIdentifier = new LeafIdentifier(nucleotide.getPdbIdentifier(), nucleotide.getModelIdentifier(), nucleotide.getChainIdentifier(), nucleotide.getSerial() - 1);
                                OakNucleotide container = ((OakNucleotide) structure.getLeafSubstructure(leafIdentifier).get());
                                container.addAtom(oxygenOne);
                                //nucleotide.addBondBetween((OakAtom) backbonePosphateOptional.get(), oxygenOne);


                            } else {
                                backboneFailCount++;
                                if (backboneFailCount == 1) {
                                    logger.warn("Could not calculate backbone for nucleotide {}.", nucleotide);
                                } else {
                                    logger.warn("Could not calculate backbone more than once.");
                                }
                            }

                            // remove obsolete atoms
                            firstPhosphateOptional.ifPresent(atom -> nucleotide.removeAtom(atom.getAtomIdentifier()));
                            secondPhosphateOptional.ifPresent(atom -> nucleotide.removeAtom(atom.getAtomIdentifier()));
                            backboneOxygenFourPrimeOptional.ifPresent(atom -> nucleotide.removeAtom(atom.getAtomIdentifier()));
                            //chain.removeNode(nucleotide.getAtomByName(AtomName.O4Pr));

                            // remove obsolete RNA specific atoms
                            backboneOxygenTwoPrimeOptional.ifPresent(atom -> nucleotide.removeAtom(atom.getAtomIdentifier()));


                            NucleotideValidator validator = new NucleotideValidator(nucleotide.getIdentifier());
                            nucleotide.getAllAtoms().forEach(atom -> convertAtom(structure,nucleotide, atom, validator));
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

    private static void convertAtom(Structure structure, OakNucleotide nucleotide, Atom an, NucleotideValidator validator) {
        switch (an.getAtomName()) {
            case "O5'": {
                String replacement = "N1'";
                replace(nucleotide, an, NITROGEN, replacement);
                validator.validate(replacement);
                break;
            }
            case "C5'": {
                String replacement = "C2'";
                replace(nucleotide, an, CARBON, replacement);
                validator.validate(replacement);
                break;
            }
            case "C4'": {
                String replacement = "C3'";
                replace(nucleotide, an, CARBON, replacement);
                validator.validate(replacement);
                break;
            }
            case "C3'": {
                String replacement = "N4'";
                replace(nucleotide, an, NITROGEN, replacement);
                validator.validate(replacement);
                break;
            }
            case "C2'": {
                String replacement = "C7'";
                replace(nucleotide, an, CARBON, replacement);
                validator.validate(replacement);
                break;
            }
            case "C1'": {
                String replacement = "C8'";
                replace(nucleotide, an, CARBON, replacement);
                validator.validate(replacement);
                break;
            }
            case "O3'": {
                String replacement = "C5'";
                replace(nucleotide, an, CARBON, replacement);
                validator.validate(replacement);
                break;
            }
            case "P": {
                String replacement = "C'";
                an = replace(nucleotide, an, CARBON, replacement);
                nucleotide.removeAtom(an.getAtomIdentifier());
                LeafIdentifier leafIdentifier = new LeafIdentifier(nucleotide.getPdbIdentifier(), nucleotide.getModelIdentifier(), nucleotide.getChainIdentifier(), nucleotide.getSerial() - 1);
                OakNucleotide container = ((OakNucleotide) structure.getLeafSubstructure(leafIdentifier).get());
                container.addAtom((OakAtom) an);

                validator.validate(replacement);
                break;
            }
            default:
                break;
        }
    }

    private static Atom replace(OakNucleotide nucleotide, Atom atom, Element replacementElement, String replacementString) {
        logger.trace("Replacing Atom {} through {}.", atom.getAtomName(), replacementString);
        OakAtom replacementAtom = new OakAtom(atom.getAtomIdentifier(), replacementElement, replacementString, atom.getPosition());
        nucleotide.removeAtom(atom.getAtomIdentifier());
        nucleotide.addAtom(replacementAtom);
        return replacementAtom;
    }

    /**
     * @param a       first Atom for the calculation of the centroid between Atom a and Atom b
     * @param b       second Atom for the calculation of the centroid between Atom a and Atom b
     * @param c       an Atom that is bounded with the calculated missing Atom
     * @param forward an boolean that defined the direction of the new normalized vector (true is positive and false
     *                negate the vector)
     * @param name    an String of the missing atom
     * @param element an Element equals the missing atom
     */
    public static OakAtom calculateMissingAtoms(Atom a, Atom b, Atom c, Boolean forward, String name, Element element) {

        Vector3D positionA = a.getPosition();
        Vector3D positionB = b.getPosition();

        Vector3D centroid = Vectors3D.getCentroid(Arrays.asList(positionA, positionB));
        logger.trace("Calculated centroid between {} and {} as {}.", a, b, centroid);

        Vector3D missingAtomPosition;
        if (forward) {
            missingAtomPosition = c.getPosition().add(centroid.subtract(c.getPosition()).normalize().multiply(C_O_DOUBLE_BOND_DISTANCE));

        } else {
            missingAtomPosition = c.getPosition().add(centroid.subtract(c.getPosition()).normalize().additivelyInvert().multiply(C_O_DOUBLE_BOND_DISTANCE));

        }

        return new OakAtom(addedAtomIndex++, element, name, missingAtomPosition);
    }
}
