# Appointment-Scheduling-Problem

It is common practice for a multitude of service providers, especially in healthcare, to
set up appointment systems for their activities. Clearly, there are numerous reasons to
do so; considerable improvements in efficiency and cost saving can be achieved. For instance, appointments allow service providers with the possibility to allocate customers
in future time slots, thereby preventing them from missing out on any job opportunity. In the context of surgery scheduling, efficient appointment scheduling leads to
improved resource utilization and consequent cost saving. On the customer’s side,
appointments are meant to considerably reduce waiting time, which is a crucial factor
that determines the quality of service. Although it is generally easy to devise basic appointment systems, obtaining optimal results often presents more than a challenge. In
particular, when uncertainty is involved, ideal outcomes cannot be achieved anymore.
For instance, it frequently happens that customers experience waiting time regardless
of a well scheduled appointment scheme. In order to work with uncertainty, concepts
from stochastic programming have been employed for constructing the mathematical
model presented in the following sections. When dealing with appointment systems,
it is important to make a distinction between scheduling and sequencing. Sequencing
refers to the branch of optimization that determines the optimal order in which customers need to be scheduled. Scheduling focuses on establishing the optimal length
of the appointments. In this project we propose a solution to the so-called 'appointment scheduling problem' with particular emphasis on the method known as ’sample
average approximation’.

We formulate a basic model and an extended model that accounts for no-shows. In either case, we developed solutions for discrete and continuous distributions.
Refer to 'ASP.pdf' for further explanations. The problems were solved by brute-force and by means of a decomposition technique. This repository contains the MATLAB code
for the two approaches applied to the basic and the extended model.
