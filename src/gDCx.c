/*
 ============================================================================
 Name        : gDCx.c
 Author      : Peter Jensen
 Version     : 0.0.1
 Copyright   : Copyright (c) 2010 Peter Jensen All Rights Reserved
 Description : Achieves 75% compression for genetic data where bytes belong
 to the set { T, C, G, A } by using a permutation table.
 ============================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

// The size of each permutation, do not change.
#define PERMUTATION_SIZE			4
#define OUTPUT_EXTENSION			".cx"

// Used for constructing permutation table
enum amino_acid { T = 0x00, C, G, A};

// Supports amino acid permutations and lookup by index.
char aminoAcidTable[256][PERMUTATION_SIZE];

// Common location for storing specific permutations.
char acid[PERMUTATION_SIZE] = {'T', 'C', 'G', 'A'};

// Structure for storing input parameters
struct parameters
{
	char fileName[25];
} params;

// Packed and unpacked compression structure.
struct compressionFooter
{
	// The total uncompressed size.
	size_t uncompressedSize, compressedSize;

	// Everything that doesn't meet the MOD 4 requirement.
	unsigned char overflow[3];
};

// Declare a place for genetic data to be stored.
unsigned char* geneticData;

// Location for compressed genetic data.
unsigned char* compressedGeneticData;

// Prototypes
void showHelp(void);

/**
 * Function: doError
 * Description: Prints an error message to standard error.
 * Parameters: aErrorMessage - The error message to put to standard error.
 */
void doError(const char* aErrorMessage)
{
	fprintf(stderr, "%s", aErrorMessage);
	exit(EXIT_FAILURE);
}

/**
 * Function: getPermutation
 * Parameters:
 * 	index - The index of the permutation to set.
 * 	dest - The destination to write the permutation to.
 * Preconditions: Assumes that the permutation table has already been
 * constructed and exists.
 */
void getPermutation(unsigned short int index, char* dest)
{
	unsigned char tmpCopy[5] = {'\0','\0','\0','\0','\0'};

	unsigned short int tmpIndex = 0x00;

	for (tmpIndex = 0x00; tmpIndex < (PERMUTATION_SIZE); ++tmpIndex)
		tmpCopy[tmpIndex] = aminoAcidTable[index][tmpIndex];

	dest[0] = '\0';

	strcat((char*)dest, (char*)tmpCopy);
}

/**
 * Function: getPermuation
 * Inputs
 * 	aPermStr - The permutation string to look for.
 * Returns: The index to that permutation in the permutation table.
 */
short int getPermuationIndex(const char* aPermStr)
{
	// Return -1 when input length does not equal permute size.
	if (strlen((char*)aPermStr) != PERMUTATION_SIZE)
		return -1;

	// Index through the values in the permutation table
	// looking for the passed permutation string.
	unsigned short int index = 0;
	for (index = 0; index <= UCHAR_MAX; ++index)
		if (aminoAcidTable[index][0] == aPermStr[0] &&
			aminoAcidTable[index][1] == aPermStr[1] &&
			aminoAcidTable[index][2] == aPermStr[2] &&
			aminoAcidTable[index][3] == aPermStr[3])
			return index;

	return -1;
}

/**
 * Function: compressGeneticData
 * Compresses the contents of geneticData by performing a lookup using the
 * permutation table and writing that index to an outfile for later
 * decompression.  Lastly, anything that doesn't match MOD 4 requirement
 * simply gets appended to the file in the form of a packed struct.
 * Inputs: aDataSize - The length of the data to compress.
 * Returns: Void
 * Preconditions: geneticData exists, permutation table exists
 * Postconditions: Compressed data written to output file.
 */
void compressGeneticData(const long aDataSize)
{
	// Points to the eventual output file.
	FILE* fp;

	// Set up two integers used strictly for indexing.
	register unsigned short int subIndex = 0;
	register unsigned short int index = 0;

	// The index of the permutation lookup.
	short int compressionIndex = 0;

	// this compression structure will be packed into fp using fwrite.
	struct compressionFooter footer;

	// The output fileName.
	char outFileName[25];

	// Show that this operation has started.
	printf("Compressing ... ");

	fflush(stdout);

	// Generate the output file name by basing it off of the input file name.
	strcpy(outFileName, params.fileName);
	strcat(outFileName, OUTPUT_EXTENSION);

	// Open the output file
	if (!(fp = fopen(outFileName, "w")))
 		doError("Failed opening out-file.");

	// Set what is the uncompressed size of the file in the header and
	// initialize what will be the compressed size to zero.
	footer.uncompressedSize = (size_t) aDataSize;
	footer.compressedSize = 0x0000;

	for (index = 1; index <= aDataSize; ++index)
	{
		if (index % PERMUTATION_SIZE == 0)
		{
			// Assign the last index of the acid array
			acid[subIndex] = geneticData[index - 1];

			// Try to get the permutation index;
			// -1 means it wasn't found possibly a bad character or newline.
			if ((compressionIndex = getPermuationIndex(acid)) != -1) {
				if (compressionIndex >= 0 && compressionIndex < 256) {
					fwrite(&compressionIndex, sizeof(short int), 1, fp);
					footer.compressedSize += 1;
				}
				else doError("Genetic permutation not found in lookup table");
			}

			subIndex = 0;
			acid[4] = '\0';
		} else { acid[subIndex++] = geneticData[index - 1]; }
	}

	// populate overflow array with what is left over
	index = 0;
	while (index++ < subIndex)
		footer.overflow[index] = *(acid + index);

	// Pack the footer structure
	fwrite(&footer,sizeof(footer),1,fp);

	fclose(fp);

	// Completed
	printf("[Done]\n");

	// Show the output file name
	printf("Output file: %s\n", outFileName);
}

/**
 * Function: buildTablePermute
 * Description:  Finds all permutations for the four basic amino acids
 * and sets them in the aminoAcidTable.
 * Affects: Globally scoped aminoAcidTable.
 * Inputs: void
 * Returns: void
 * Preconditions: aminoAcidTable has been properly allocated.
 * Postconditions: aminoAcidTable has been fully populated.
 */
void buildTablePermute()
{
	// char acid[4] = {'T', 'C', 'G', 'A'};

	unsigned char a, b, c, d;
	unsigned short int index = 0;

	printf("Building permutation table ... ");

	// Building the permutation table
	for (a = T; a <= A; ++a)
		for (b = T; b <= A; ++b)
			for (c = T; c <= A; ++c)
				for (d = T; d <= A; ++d, ++index)
				{
					// Assign amino acids for this inner iteration.
					aminoAcidTable[index][0] = acid[a];
					aminoAcidTable[index][1] = acid[b];
					aminoAcidTable[index][2] = acid[c];
					aminoAcidTable[index][3] = acid[d];

					// Increment index embedded in last for loop
				}
	printf("[Done]\n");
}

/**
 * Function: readGeneticData
 * Description: Reads genetic data
 * Returns: The number of bytes in the data sequence.
 * Postconditions:
 *	The value of geneticData will be populated with the genetic data.
 */
const long readGeneticData()
{
	FILE* fp;				// Input file pointer
	long tmpFileSize;		// The size of the file; figure out before hand.
	size_t tmpReadResult;	// The result of doing fread on the file.
	int index;

	// Show that were about to read in genetic data.
	printf("Reading genetic data ... ");

	// Try to open params.fileName, exit on failure.
	if ((fp = fopen(params.fileName, "r")) == NULL)
		doError("Input file does not exist.");

	// Determine the file's size, set as a long value.
	fseek(fp, 0, SEEK_END);
	tmpFileSize = ftell(fp);
	rewind(fp);

	// Allocate a chunk of memory that is exactly the size of the input file.
	geneticData = (unsigned char*) malloc(sizeof(unsigned char) * tmpFileSize);

	//  Make sure we didn't run out of memory.
	if (geneticData == NULL)
		doError("Out of memory.");

	// Read genetic data from input file.
	tmpReadResult = fread(geneticData, 1, tmpFileSize, fp);

	// make sure the size of what was read and expected values match.
	if (tmpReadResult != tmpFileSize)
		doError("Input file read failure.");

	fclose(fp);

	// Make sure the data ends up as upper case.
	for (index = 0; index < tmpFileSize; ++index)
		*(geneticData + index) = toupper(geneticData[index]);

	printf("[Done]\n");

	return tmpFileSize;
}

/**
 * Function: parseArgs
 * Description: Parse arguments and populate the parameters struct with
 * settings pulled from argv.
 * Remarks: argv and argc are not modified.
 */
void parseArgs(const int argc, char** argv)
{
	// Define various parameters for command line usage.
	const char* FILENAME_PARAM = "-f";
	const char* HELP_PARAM = "-help";

	int index;

	// Index through all arguments
	for (index = 0; index < argc; ++index)
		if (strcmp(FILENAME_PARAM, argv[index]) == 0) {
			if (argv[index + 1] != '\0')
				strcat(params.fileName, argv[index + 1]);
		} else if (strcmp(HELP_PARAM, argv[index]) == 0) {
			(void) showHelp();
			exit(EXIT_SUCCESS);
		}

	// Show various settings.
	printf("\n\nFilename: %s\n", params.fileName);
}

void showHelp() {
	printf("gDCx :: Genetic Data Compression Utility");
	printf("<Usage>");
	printf("   ./gDCx -f <filename> -help");
}


int main(int argc, char** argv) {

	// Stores the size of the genetic code.
	long tmpSize;

	// Parses incoming arguments like fileName.
	parseArgs(argc, argv);

	// Read genetic data will return the long size of the data.
	tmpSize = readGeneticData();

	// Build the permutation table.
	buildTablePermute();

	// Compress and write genetic data
	compressGeneticData(tmpSize);

	free(geneticData);

	return EXIT_SUCCESS;
}




