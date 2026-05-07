# MHGuide UMR

Parser-Library zum Einlesen von MHGuide JSON Dateien unter Berücksichtigung der Besonderheiten am Standort Marburg.

## RNA Fusionen

RNA Fusionen werden exportiert, wenn die Angaben in der JSON-Datei unter `REPORT_NARRATIVE` vorhanden sind.
Hierbei ist das Format einzuhalten:

* Jede Fusion wird in einer neuen Zeile gelistet.
* Jede Teilangabe der Fusion wird durch ein Semikolon getrennt.

Beispiel:
`ABCD1(ex 1)::ABCD2(ex 2), transcript ID: NM_012345.4/NM_012456.2, strand: -/-, breakpoint: chr19:12345678/chr19:13456789, supporting read pairs: 1234`
