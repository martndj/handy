Human And Nature DYnamical model
================================

HANDY model was formulated in  
Motesharrei S., Rivas, J. and Kalnay E., __Human and Nature Dynamics (HANDY): Modeling Inequality and Use of Resources in the Collapse or Sustainability of Societies__, 2014

Here lies my implementation of HANDY!  
_May we evolve before the collapse!_

Martin Deshaies-Jacques

[deshaies.martin@sca.uqam.ca](mailto:deshaies.martin@sca.uqam.ca)

[www.science.martn.info](http://www.science.martn.info)

Abstract from Motesharrei et al. (2014)
---------------------------------------
There are widespread concerns that current trends in resource-use are unsustainable, but possibilities of overshoot/collapse remain controversial. Collapses have occurred frequently in history, often followed by centuries of economic, intellectual, and population decline. Many different natural and social phenomena have been invoked to explain specific collapses, but a general explanation remains elusive.

In this paper, we build a human population dynamics model by adding accumulated wealth and economic inequality to a predator-prey model of humans and nature. The model structure, and simulated scenarios that offer significant implications, are explained. Four equations describe the evolution of Elites, Commoners, Nature, and Wealth. The model shows Economic Stratification or Ecological Strain can independently lead to collapse, in agreement with the historical record.

The measure “Carrying Capacity” is developed and its estimation is shown to be a practical means for early detection of a collapse. Mechanisms leading to two types of collapses are discussed. The new dynamics of this model can also reproduce the irreversible collapses found in history. Collapse can be avoided, and population can reach a steady state at maximum carrying capacity if the rate of depletion of nature is reduced to a sustainable level and if
resources are distributed equitably.

Use case
--------
Put handy folder in your $PYTHONPATH.

    import handy
    m=handy.handy(depletionR=20.e-6)
    m.integrate(900, x0=[100., 25.])
    ax1,ax2=handy.plotState(m)
    ax1.legend(loc='upper left')
    ax2.legend(loc='upper right')
