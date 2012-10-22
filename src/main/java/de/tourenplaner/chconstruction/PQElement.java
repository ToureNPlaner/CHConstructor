/*
 * (C) Copyright 2012 Peter Vollmer
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; version 2
 * of the License.
 */

package de.tourenplaner.chconstruction;

class PQElement implements Comparable<PQElement> {

    public int key;
    public int value;

    PQElement(int a, int b) {
        key = a;
        value = b;
    }

    public int compareTo(PQElement o) {
        if (key > o.key) return 1;
        else if (key == o.key) return 0;
        else return -1;
    }
}