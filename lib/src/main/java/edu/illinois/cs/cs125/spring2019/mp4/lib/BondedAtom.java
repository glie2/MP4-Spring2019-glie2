package edu.illinois.cs.cs125.spring2019.mp4.lib;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;

/**
 * Represents an atom in an organic molecule and its bonds.
 * <p>
 * Since references to all other atoms in the molecule can be reached through the methods of this class, it also
 * represents the overall molecule.
 * <p>
 * One atom can be connected to another through at least one bond. For example, carbon should have 4 bonds, but may
 * be connected to only 3 other atoms if one of the connections is a double bond.
 * <p>
 * For convenience BondedAtom implements the Iterable interface allow you to iterate through all of its neighbors
 * using the enhanced Java for-each loop.
 */
public final class BondedAtom implements Iterable<BondedAtom> {
    /**
     * What element this atom is of.
     */
    private ChemicalElement element;
    /**
     * The atoms that this one is connected to.
     */
    private List<BondInfo> connectedAtoms;
    /**
     * Reference to an array of bonds created during molecule construction.
     */
    private BondedAtom[] otherAtoms;
    /**
     * Reference to an array of bound counts created during molecule construction.
     */
    private int[] bondCounts;

    /**
     * The ID of this instance, for ease of debugging.
     */
    private int id;
    /**
     * The last assigned instance ID.
     */
    private static int lastId = 0;

    /**
     * Creates a new atom connected to the given others.
     *
     * @param setElement    What element this atom is.
     * @param setOtherAtoms The distinct other atoms this atom is connected to.
     * @param setBondCounts The bond counts in each connection.
     */
    public BondedAtom(final ChemicalElement setElement, final BondedAtom[] setOtherAtoms,
                      final int[] setBondCounts) {
        id = ++lastId;
        element = setElement;

        if (setOtherAtoms == null) {
            return;
        }
        /*
         * These don't make copies because otherwise it would be impossible to create an atom bonded
         * to others, since those would have to be created first, but creating those would require
         * that this is already created!
         */
        otherAtoms = setOtherAtoms;
        bondCounts = setBondCounts;
        /*
         * Unfortunately we can't immediately create our array of BondInfo objects here, since
         * the arrays that we are passed are still being filled in by the molecule creator.
         * Instead we'll defer this until the first time that the BondedAtom is read.
         */
    }

    /**
     * Gets what element this atom is an instance of.
     *
     * @return This atom's element.
     */
    public ChemicalElement getElement() {
        return element;
    }

    /**
     * Helper function to load our array of BondInfo objects.
     */
    private void loadConnectedAtoms() {
        if (connectedAtoms != null) {
            return;
        }
        connectedAtoms = new ArrayList<>(otherAtoms.length);
        for (int i = 0; i < otherAtoms.length; i++) {
            connectedAtoms.add(new BondInfo(otherAtoms[i], bondCounts[i]));
        }
    }

    /**
     * Return an iterator of BondedAtoms connected to this BondedAtom.
     * <p>
     * This is commonly used by the molecular analysis code to explore the molecule graph.
     *
     * @return an iterator of BondedAtoms connected to this BondedAtom.
     */
    public Iterator<BondedAtom> iterator() {
        loadConnectedAtoms();
        List<BondedAtom> atoms = new ArrayList<>();
        for (BondInfo bond : connectedAtoms) {
            atoms.add(bond.getAtom());
        }
        return atoms.iterator();
    }

    /**
     * Return an iterable list of bond information about all neighboring atoms.
     * <p>
     * Each BondInfo object includes both a reference to the neighbor and it bond count.
     * It is frequently more convenient to use the other iterator method above that returns BondedAtoms.
     * Only use this if you need <i>both</i> the BondedAtom <i>and</i> the bond count.
     *
     * @return an unmodifiable list of BondInfo objects describing our neighbor atoms.
     */
    public List<BondInfo> getBondInfo() {
        loadConnectedAtoms();
        return Collections.unmodifiableList(connectedAtoms);
    }

    /**
     * Convenience method used by the test suite.
     *
     * @param index index of the BondedAtom to return
     * @return the BondedAtom at that position or null if the index is invalid
     */
    public BondedAtom getConnectedAtom(final int index) {
        loadConnectedAtoms();
        try {
            return connectedAtoms.get(index).getAtom();
        } catch (IndexOutOfBoundsException unused) {
            return null;
        }
    }

    /**
     * Convenience method for determining if this atom is carbon.
     *
     * @return whether this atom is carbon
     */
    public boolean isCarbon() {
        return element == ChemicalElement.CARBON;
    }

    /**
     * Convenience method for determining if this atom is a substituent.
     *
     * @param backbone the backbone of this molecule to use to identify off-backbone carbons
     * @return whether this atom is a substituent
     */
    public boolean isSubstituent(final List<BondedAtom> backbone) {
        return (element != ChemicalElement.CARBON && element != ChemicalElement.HYDROGEN)
            || (element == ChemicalElement.CARBON && !(backbone.contains(this)));
    }

    /**
     * Convenience method for determining of this atom has a substituent.
     * <p>
     * This could either mean it has an attached substituent or off-backbone carbon.
     *
     * @param backbone the backbone of this molecule to use to identify off-backbone carbons
     * @return true if this atom has a connected substituent, and false otherwise
     */
    public boolean hasSubstituent(final List<BondedAtom> backbone) {
        loadConnectedAtoms();

        for (BondedAtom connected : this) {
            if (connected.isSubstituent(backbone)) {
                return true;
            }
        }
        return false;
    }

    /**
     * Return the count of high-priority substituents attached to this atom.
     *
     * @return the count of high-priority substituents attached to this atom
     */
    public int highPrioritySubstituentCount() {
        loadConnectedAtoms();

        int oxygens = 0;
        for (BondedAtom connected : this) {
            if (connected.getElement() == ChemicalElement.OXYGEN) {
                // This double-counts carboxylic acids (=O and -OH)
                // But that doesn't matter because they're going to be carbon #1 no matter what
                oxygens++;
            }
        }
        return oxygens;
    }

    /**
     * Return the count of low-priority substituents attached to this atom.
     *
     * @param backbone the backbone of this molecule to use to identify off-backbone carbons
     * @return the count of low-priority substituents attached to this atom
     */
    public int lowPrioritySubstituentCount(final List<BondedAtom> backbone) {
        loadConnectedAtoms();

        int substituents = 0;
        for (BondedAtom connected : this) {
            if (connected.getElement() == ChemicalElement.CARBON && !(backbone.contains(connected))) {
                substituents++;
            } else if (NamingConstants.HALIDES.containsKey(connected.getElement())) {
                substituents++;
            }
        }
        return substituents;
    }

    /**
     * Gets the name of a substituent starting at this carbon.
     * <p>
     * Assumes that there is no end group on the core/base carbon.
     *
     * @param substituent the first atom of the substituent
     * @return the substituent's name, or null if none exists or this atom is not carbon
     */
    public String nameSubstituent(final BondedAtom substituent) {
        if (!(isCarbon())) {
            return null;
        }

        loadConnectedAtoms();

        if (substituent.getElement() == ChemicalElement.CARBON) {
            List<BondedAtom> path = new ArrayList<>(Arrays.asList(substituent, this));
            path:
            while (true) {
                BondedAtom previous = path.get(1);
                for (BondedAtom connected : path.get(0)) {
                    if (connected != previous && connected.getElement() == ChemicalElement.CARBON) {
                        path.add(0, connected);
                        continue path;
                    }
                }
                break;
            }
            return NamingConstants.CHAIN_BASE_NAMES[path.size() - 2] + "yl";
        } else if (substituent.getElement() == ChemicalElement.OXYGEN) {
            for (BondedAtom.BondInfo bondInfo : this.getBondInfo()) {
                if (bondInfo.getAtom() != substituent) {
                    continue;
                }
                if (bondInfo.getCount() == 2) {
                    return "one";
                } else {
                    return "ol";
                }
            }
            return null;
        } else {
            return NamingConstants.HALIDES.get(substituent.getElement());
        }
    }

    /**
     * Gets the name of the end group substituent, if any, on this carbon.
     *
     * @return the name of any end substituent, or null if none exists or this atom is not carbon
     */
    public String nameEndGroup() {
        loadConnectedAtoms();

        if (!(isCarbon())) {
            return null;
        }
        boolean hasCarbonyl = false;
        boolean hasHydroxyl = false;
        for (BondedAtom.BondInfo bondInfo : this.getBondInfo()) {
            if (bondInfo.getAtom().getElement() == ChemicalElement.OXYGEN) {
                if (bondInfo.getCount() == 2) {
                    hasCarbonyl = true;
                } else {
                    hasHydroxyl = true;
                }
            }
        }
        if (hasCarbonyl) {
            if (hasHydroxyl) {
                return "oic acid";
            } else {
                return "al";
            }
        } else {
            return null;
        }
    }

    /**
     * Returns the string representation of this atom.
     *
     * @return This atom's element and instance ID.
     */
    @Override
    public String toString() {
        return getElement().getSymbol() + " #" + id;
    }

    /**
     * Class storing information about a bond to one of our connectedAtoms.
     * <p>
     * Stores both a reference to the neighbor itself and the count of bonds to that neighbor.
     */
    class BondInfo {
        /**
         * Neighbor this atom is bonded to.
         */
        private BondedAtom atom;

        /**
         * Return the neighbor this atom is bonded to.
         *
         * @return the neighbor this atom is bonded to.
         */
        public BondedAtom getAtom() {
            return atom;
        }

        /**
         * Count of bonds to that neighbor.
         */
        private int count;

        /**
         * Return the count of bonds to this neighbor.
         *
         * @return the count of bonds to this neighbor.
         */
        public int getCount() {
            return count;
        }

        /**
         * Create a new BondInfo object to store data about a bond to one of our connectedAtoms.
         *
         * @param setAtom  neighbor this atom is bonded to
         * @param setCount count of bounds to this neighbor
         */
        BondInfo(final BondedAtom setAtom, final int setCount) {
            atom = setAtom;
            count = setCount;
        }

        /**
         * Returns the string representation of this BondInfo object.
         *
         * @return a string indicating the connected atom and the bond strength
         */
        @Override
        public String toString() {
            return "{BondInfo " + "-=%".substring(count - 1, count) + atom.getElement().getSymbol() + "}";
        }
    }
}
