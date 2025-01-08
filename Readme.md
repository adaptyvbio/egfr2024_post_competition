# EGFR2024 post competition data & script collection repo

Repo organized by name of team/org/contributor, e.g. `./adaptyv/` will contain scripts and data added by adaptyv, with `./adaptyv/all_submissions.parquet` being a dataframe of all submission data.

# Git LFS Setup and Usage Guide

This repository uses Git Large File Storage (LFS) to handle large "binary" files ala csv, pdf,parquet etc.

## Setting Up Git LFS

1. Install Git LFS:

   **For Mac (using Homebrew):**

   ```bash
   brew install git-lfs
   ```

   **For Windows:**

   - Download and install from https://git-lfs.github.com

   **For Linux (Debian/Ubuntu):**

   ```bash
   sudo apt-get install git-lfs
   ```

2. Initialize Git LFS:
   ```bash
   git lfs install
   ```

## Cloning the Repository

1. Clone the repository as normal:

   ```bash
   git clone git@github.com:adaptyvbio/egfr2024_post_competition.git
   ```

2. The LFS files will be downloaded automatically during the clone process.

## Working with LFS Files

- When you add new files of the tracked types ( e.g. _.csv, _.pdf, \*.parquet), they will automatically be handled by Git LFS.
- You can verify if a file is being tracked by LFS using:
  ```bash
  git lfs ls-files
  ```

## Pulling LFS Files

If you need to pull the latest LFS files:

```bash
git lfs pull
```

## Checking LFS Status

To see the status of your LFS files:

```bash
git lfs status
```

## Troubleshooting

If you encounter issues:

1. Verify Git LFS is installed:

   ```bash
   git lfs version
   ```

2. Ensure LFS is initialized in your repository:

   ```bash
   git lfs install
   ```

3. If files aren't being tracked properly, verify the .gitattributes file contains:
   ```
   *.csv filter=lfs diff=lfs merge=lfs -text
   *.pdf filter=lfs diff=lfs merge=lfs -text
   *.parquet filter=lfs diff=lfs merge=lfs -text
   ```
   etc.

## Additional Commands

- Track a new file pattern:

  ```bash
  git lfs track "*.extension"
  ```

- Untrack a file pattern:
  ```bash
  git lfs untrack "*.extension"
  ```

For more information, visit the [Git LFS website](https://git-lfs.github.com/).

```

This README provides comprehensive instructions for:
1. Installing Git LFS on different operating systems
2. Setting up Git LFS in a repository
3. Working with LFS files
4. Troubleshooting common issues
5. Additional useful commands

Users can follow these instructions to properly work with the large files in your repository.
```
