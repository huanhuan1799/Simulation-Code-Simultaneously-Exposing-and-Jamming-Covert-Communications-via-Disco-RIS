# Simulation Code: Simultaneously Exposing and Jamming Covert Communications via Disco RIS

This repository contains the MATLAB simulation code for the paper: **"Simultaneously Exposing and Jamming Covert Communications via Disco Reconfigurable Intelligent Surfaces"**, published in **IEEE Journal on Selected Areas in Communications (JSAC)**.

We highly respect reproducible research and aim to provide open access to the simulation codes used in our published work.

## 📄 Abstract

Covert communications provide a stronger privacy protection than cryptography and physical-layer security (PLS). However, previous works on covert communications have implicitly assumed the validity of channel reciprocity, i.e., wireless channels remain constant or approximately constant during their coherence time. In this work, we investigate covert communications in the presence of a disco RIS (DRIS) deployed by the warden Willie, where the DRIS with random and time-varying reflective coefficients acts as a “disco ball”, introducing time-varying fully-passive jamming (FPJ). Consequently, the channel reciprocity assumption no longer holds. The DRIS not only jams the covert transmissions between Alice and Bob, but also decreases the error probabilities of Willie’s detections, without either Bob’s channel knowledge or additional jamming power. To quantify the impact of the DRIS on covert communications, we first design a detection rule for the warden Willie in the presence of time-varying FPJ introduced by the DRIS. Then, we define the detection error probabilities, i.e., the false alarm rate (FAR) and the missed detection rate (MDR), as the monitoring performance metrics for Willie’s detections, and the signal-to-jamming-plus-noise ratio (SJNR) as a communication performance metric for the covert transmissions between Alice and Bob. Based on the detection rule, we derive the detection threshold for the warden Willie to detect whether communications between Alice and Bob is ongoing, considering the time-varying DRIS-based FPJ. Moreover, we conduct theoretical analyses of the FAR and the MDR at the warden Willie, as well as SJNR at Bob, and then present unique properties of the DRIS-based FPJ in covert communications. We present numerical results to validate the derived theoretical analyses and evaluate the impact of DRIS on covert communications.

## 🚀 Usage

To reproduce the simulation results presented in the paper (specifically **Fig. 3**), please follow these steps:

1.  Ensure you have MATLAB installed.
2.  Download or clone this repository.
3.  Run the main script:
    ```matlab
    JSAC_CovertCommun.m
    ```

**Note:** The file `Covert_Communication_GenerateAJchannel.m` is a helper function used to generate the wireless channel realizations. All simulation results presented in the paper can be obtained by slightly modifying the parameters in these two files.

## 📂 File Description

* **`JSAC_CovertCommun.m`**: The main simulation script. It handles the Monte Carlo simulations for False Alarm Rate (FAR), Missed Detection Rate (MDR), and Signal-to-Jamming-plus-Noise Ratio (SJNR).
* **`Covert_Communication_GenerateAJchannel.m`**: A function to generate the channel coefficients, including the direct links and the cascaded RIS links.

## 📖 Citation

If you use this simulation code or find our work helpful in your research, please kindly cite our paper:

**Plain Text:**
> H. Huang, H. Zhang, Y. Cai, D. Niyato, A. L. Swindlehurst, and Z. Han, "Simultaneously Exposing and Jamming Covert Communications via Disco Reconfigurable Intelligent Surfaces," *IEEE Journal on Selected Areas in Communications*, Dec. 2025. DOI: 10.1109/JSAC.2025.3646955.

**BibTeX:**
```bibtex
@article{Huang2025JSAC,
  author    = {Huang, Huan and Zhang, Hongliang and Cai, Yi and Niyato, Dusit and Swindlehurst, A. Lee and Han, Zhu},
  title     = {Simultaneously Exposing and Jamming Covert Communications via Disco Reconfigurable Intelligent Surfaces},
  journal   = {IEEE Journal on Selected Areas in Communications},
  year      = {2025},
  month     = {Dec},
  doi       = {10.1109/JSAC.2025.3646955}
}


## 📚 Related Work (The "DISCO" Series)

For more details about the **DISCO** concept and related ideas, please refer to our other published papers:

H. Huang, L. Dai, H. Zhang, C. Zhang, Z. Tian, Y. Cai, A. L. Swindlehurst, and Z. Han, “DISCO might not be funky: Random intelligent reflective surface configurations that attack,” *IEEE Wireless Communications*, vol. 31, no. 5, pp. 76-82, Oct. 2024.

H. Huang, Y. Zhang, H. Zhang, C. Zhang, and Z. Han, “Illegal intelligent reflecting surface based active channel aging: When jammer can attack without power and CSI,” *IEEE Transactions on Vehicular Technology*, vol. 72, no. 8, pp. 11018-11022, Aug. 2023.

H. Huang, Y. Zhang, H. Zhang, Y. Cai, A. L. Swindlehurst, and Z. Han, “Disco intelligent reflecting surfaces: Active channel aging for fully-passive jamming attacks,” *IEEE Transactions on Wireless Communications*, vol. 23, no. 1, pp. 806–819, Jan. 2024.

H. Huang, L. Dai, H. Zhang, Z. Tian, Y. Cai, C. Zhang, A. L. Swindlehurst, and Z. Han, “Anti-jamming precoding for disco intelligent reflecting surfaces based fully-passive jamming attacks,” *IEEE Transactions on Wireless Communications*, vol. 23, no. 8, pp. 9315-9329, Aug. 2024.

H. Huang, H. Zhang, Y. Cai, A. L. Swindlehurst, and Z. Han, “Disco Intelligent Omni-Surfaces: 360° Fully-Passive Jamming Attacks,” *IEEE Transactions on Wireless Communications*, early access, Jun. 2025. DOI: 10.1109/TWC.2025.3581208.
