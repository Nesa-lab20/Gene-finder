
#main function 
import argparse

def main():
    parser = argparse.ArgumentParser(description="Gene Finder Tool")
    parser.add_argument("input_file", type=str, help="Path to the FASTA file")
    args = parser.parse_args()

    # Call function to find genes (implementation pending)
    print(f"Reading file: {args.input_file}")

if __name__ == "__main__":
    main()
