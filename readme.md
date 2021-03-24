{
  "nbformat": 4,
  "nbformat_minor": 0,
  "metadata": {
    "colab": {
      "name": "readme_aa_grouping.ipynb",
      "provenance": [],
      "authorship_tag": "ABX9TyN7Ls583jlqbbWyawFFV5pH",
      "include_colab_link": true
    },
    "kernelspec": {
      "name": "python3",
      "display_name": "Python 3"
    },
    "language_info": {
      "name": "python"
    }
  },
  "cells": [
    {
      "cell_type": "markdown",
      "metadata": {
        "id": "view-in-github",
        "colab_type": "text"
      },
      "source": [
        "<a href=\"https://colab.research.google.com/github/TitliSarkar/Amino-Acid-Grouping/blob/main/readme.md\" target=\"_parent\"><img src=\"https://colab.research.google.com/assets/colab-badge.svg\" alt=\"Open In Colab\"/></a>"
      ]
    },
    {
      "cell_type": "markdown",
      "metadata": {
        "id": "dHXVhrdMVRQg"
      },
      "source": [
        "**Amino-Acid-Grouping**\n",
        "\n",
        "This repository contains the data and codes for TSR key representation of pdbs with and without amino acid grouping. This repositary also inclues codes for pdb sample generastion, pairwise global similarity calculation (Jaccard) and local similarity calculation (patch-based method) between proteins using TSR keys, generation of common keys and visualization.\n"
      ]
    },
    {
      "cell_type": "markdown",
      "metadata": {
        "id": "hOij-C8zVh2Y"
      },
      "source": [
        "**Description of the folder structure:**\n",
        "\n",
        "*/data* - Data directory contains all the dataset we have used in experements and details of a sample dataset.\n",
        "\n",
        "*/codes* - The codes are listed in the order of how they needs to be executed.\n",
        "      /*1_sample-generation* - Parallel code for retrieving pdb files from DB databank https://www.rcsb.org/.\n",
        "      */2_TSR-key-transformation *- Code for transforming pdb to TSR keys, with and without amino acid grouping."
      ]
    },
    {
      "cell_type": "markdown",
      "metadata": {
        "id": "gy0AGJAHWjfF"
      },
      "source": [
        ""
      ]
    }
  ]
}