#!/usr/bin/env bash

# Copy this script to the build folder

make -j 8 
./TestReferenceAmoebaMultipoleForce
./TestReferenceMBPolOneBodyForce
./TestReferenceMBPolTwoBodyForce
./TestReferenceMBPolThreeBodyForce
./TestReferenceMBPolDispersionForce
./TestReferenceMBPolIntegrationTest
