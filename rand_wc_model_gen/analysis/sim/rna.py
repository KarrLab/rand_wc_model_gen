""" Analysis of RNA dynamics of simulations of random whole-cell models

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2018-06-14
:Copyright: 2018, Karr Lab
:License: MIT
"""
import glob
import numpy
import os
import wc_analysis
import wc_sim.multialgorithm.run_results


class RnaSimulationAnalysis(wc_analysis.SimulationAnalysis):
    """ Analysis of RNA dynamics of simulations of random whole-cell models """

    def run(self):
        # load run results
        self._run_results = []
        for sim_results_path in os.listdir(self.sim_results_path):
            self._run_results.append(wc_sim.multialgorithm.run_results.RunResults(
                os.path.join(self.sim_results_path, sim_results_path)).run_results)

        # draw and show/save figures
        self.plot_individual_rna()
        self.plot_total_rna()

    def plot_individual_rna(self, i_sim=0):
        """ Plot the dynamics of each RNA molecule for a single simulation """

        # get data
        time = self._run_results[0]['populations'].index

        rna_ids = [st.id + '[c]' for st in self.model.species_types if st.id.startswith('rna_')]
        rna_cnts = self._run_results[i_sim]['populations'][rna_ids]

        # draw figure
        fig, axes = self.create_fig()
        axes.plot(time / 3600, rna_cnts.values)
        axes.set_xlim([time[0] / 3600, time[-1] / 3600])
        axes.set_ylim([0, rna_cnts.max().max()])
        axes.set_xlabel('Time (h)')
        axes.set_ylabel('RNA (molecules)')

        # show/save figure
        self.show_or_save_fig(fig, 'Individual RNA (simulation {}).pdf'.format(i_sim + 1))

    def plot_total_rna(self):
        """ Plot the dynamics of the total RNA count of each simulation """

        rna_ids = [st.id + '[c]' for st in self.model.species_types if st.id.startswith('rna_')]

        time = self._run_results[0]['populations'].index

        tot_rna_cnts = numpy.full((self._run_results[0]['populations'].shape[0], len(self._run_results)), numpy.nan)
        for i_sim, run_results in enumerate(self._run_results):
            tot_rna_cnts[:, i_sim] = run_results['populations'][rna_ids].sum(axis=1).values

        # draw figure
        fig, axes = self.create_fig()
        axes.plot(time / 3600, tot_rna_cnts)
        axes.set_xlim([time[0] / 3600, time[-1] / 3600])
        axes.set_ylim([0, numpy.max(tot_rna_cnts[:])])
        axes.set_xlabel('Time (h)')
        axes.set_ylabel('Total RNA (molecules)')

        # show/save figure
        self.show_or_save_fig(fig, 'Total RNA (all simulations).pdf')
