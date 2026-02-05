import argparse
import logging
import sys
from pathlib import Path

from simnet_hgt.pipeline import HGTPipeline


def main():
    # Default root path
    default_root = Path(__file__).parents[1] / "datasets/aminoacyl"

    parser = argparse.ArgumentParser(description="Run HGT detection pipeline.")
    parser.add_argument(
        "dataset_path",
        nargs="?",
        type=Path,
        default=default_root,
        help=f"Path to the dataset root directory (default: {default_root})",
    )
    parser.add_argument(
        "--min_coverage",
        type=float,
        default=0.8,
        help="Minimum coverage for alignment (default: 0.8)",
    )
    parser.add_argument(
        "--min_identity",
        type=float,
        default=0.59,
        help="Minimum identity for alignment (default: 0.59)",
    )

    args = parser.parse_args()

    # Configure logging
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    logger = logging.getLogger(__name__)

    try:
        pipeline = HGTPipeline(
            dataset_path=args.dataset_path,
            min_coverage=args.min_coverage,
            min_identity=args.min_identity,
        )
        pipeline.run()
    except Exception as e:
        logger.error(f"Pipeline failed: {e}", exc_info=True)
        sys.exit(1)


if __name__ == "__main__":
    main()
