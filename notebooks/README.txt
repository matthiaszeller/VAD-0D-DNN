

=== UBUNTU 18.04 VS 20.04

Upgrading to 20.04 produced errors when using :

from ipypublish import nb_setup
plt = nb_setup.setup_matplotlib(output=('pdf','svg'), rcparams=rcparams)


=> SOLUTION: install cm-super