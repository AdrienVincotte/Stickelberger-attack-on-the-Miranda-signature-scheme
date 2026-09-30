Here you will find the implementation of the attack on the Miranda signature scheme presented in the paper: https://arxiv.org/abs/2609.30925.

In the main file, you can modify the protocol parameters at the beginning of the script. Although the attack is polynomial, its complexity is determined by a formula of degree approximately 11 (assuming the linear algebra constant omega is 2.8). Unless you have an extremely powerful machine, it is not recommended to go beyond m=15.
The code includes the reduction to the MinRank instance, as well as the computation of possible gamma bases using Stickelberger's method. The final line verifies that the basis actually used to construct the keys is among the elements returned by the attack. If you wish to display all the returned bases, uncomment line 63.

Currently, the implementation is operational for cases where the minors provide enough equations to solve the system, without the need to derive new ones by increasing the degree (i.e., the d=2 case).
