#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Unit tests for the likelihood ratio bounds (Analysis.lrb()).

lrb() is deterministic (no randomness involved), so "baseline" here
means the solution/bounds computed by the pre-vectorization,
eta-direction-only implementation, captured once via
`git show HEAD:src/predictr.py` before lrb() was changed to also scan
in the beta direction (to close the gaps at the extreme left/right of
the contour) and to locally refine rows/columns with only one crossing.

Since the new implementation intentionally finds *more* solution
points (both directions + refinement), this test checks that the
solution cloud and the resulting bounds still agree closely with the
baseline, and that it does not shrink relative to it.
"""
import os
import sys
import unittest

import numpy as np
from scipy.spatial import ConvexHull

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from predictr import Analysis  # noqa: E402

DATA = [3968, 2957, 5352, 4270, 5936, 3608, 5060, 4884, 6768, 4132,
        4620, 5088, 3572, 4964, 5324]
SMALL = [3968, 2957, 5352, 4270, 5936]
CENS_DF = [4270, 5936]
CENS_DS = [3968, 2957, 5352, 3608, 5060, 4884, 6768, 4132, 4620, 5088, 3572, 4964, 5324]

# Captured once from the pre-vectorization implementation (bounds_type='2s', cl=0.9).
BASELINE = {
    'normal': {
        'n_pairs': 538,
        'beta_range': (3.758541, 7.119442),
        'eta_range': (4648.6623, 5546.5493),
        'bounds_lower': [781.5819233319696, 939.3209496171551, 1045.7115534850145],
        'bounds_upper': [1998.4192865376874, 2204.3417223059837, 2334.6353960241927],
    },
    'small_n': {
        'n_pairs': 534,
        'beta_range': (2.499477, 8.372118),
        'eta_range': (4060.8664, 5870.9095),
        'bounds_lower': [289.0432519193664, 380.5675393291786, 447.046224472208],
        'bounds_upper': [2281.5600823009268, 2480.834107954163, 2605.44120356499],
    },
    'censored': {
        'n_pairs': 532,
        'beta_range': (2.526369, 13.590700),
        'eta_range': (6047.2325, 12033.9696),
        'bounds_lower': [654.2090221039622, 860.9100469297238, 1008.8013957274502],
        'bounds_upper': [4017.381449913263, 4235.298484269525, 4369.1353646696825],
    },
}


class TestLRB(unittest.TestCase):
    def _check(self, key, x):
        base = BASELINE[key]

        # The bidirectional search + local refinement should never find
        # fewer solution points than the eta-direction-only baseline.
        self.assertGreaterEqual(len(x.beta_lrb), base['n_pairs'])

        # The overall parameter range covered by the contour should
        # match closely (small differences are expected because the
        # new mesh construction is not pixel-identical to the old one).
        # A somewhat looser tolerance here: the bidirectional search is
        # expected to nudge the extreme beta/eta values slightly beyond
        # what the eta-direction-only baseline found, since that is the
        # gap it was built to close.
        np.testing.assert_allclose(
            [x.beta_lrb.min(), x.beta_lrb.max()], base['beta_range'], rtol=2e-2)
        np.testing.assert_allclose(
            [x.eta_lrb.min(), x.eta_lrb.max()], base['eta_range'], rtol=2e-2)

        np.testing.assert_allclose(x.bounds_lower[:3], base['bounds_lower'], rtol=2e-2)
        np.testing.assert_allclose(x.bounds_upper[:3], base['bounds_upper'], rtol=2e-2)

        # The solution cloud should form a well-defined, non-degenerate
        # hull (sanity check for contour_plot's style='hull').
        hull = ConvexHull(np.column_stack((x.beta_lrb, x.eta_lrb)))
        self.assertGreater(hull.volume, 0)

    def test_lrb_uncensored(self):
        x = Analysis(df=DATA, bounds='lrb', bounds_type='2s', cl=0.9)
        x.mle()
        self._check('normal', x)

    def test_lrb_small_sample(self):
        x = Analysis(df=SMALL, bounds='lrb', bounds_type='2s', cl=0.9)
        x.mle()
        self._check('small_n', x)

    def test_lrb_censored(self):
        x = Analysis(df=CENS_DF, ds=CENS_DS, bounds='lrb', bounds_type='2s', cl=0.9)
        x.mle()
        self._check('censored', x)

    def test_lrb_no_flat_edge_near_tangent_point(self):
        """
        Regression test for a coordinate-quantization bug: near a
        steep/near-vertical stretch of the contour (here, the left tip
        at high cl), many crossings used to snap to the exact same
        coarse mesh coordinate, rendering as a spurious straight line
        segment in contour_plot's convex hull instead of the true
        curve. bisect_refine() in lrb() fixes this by refining every
        crossing to near machine precision within its bracket, so no
        coordinate value should be shared by more than a couple of
        points, and none of those should span a wide eta range.
        """
        failures_a = [0.30481336314657737, 0.5793918872111126, 0.633217732127894,
                      0.7576700925659532, 0.8394342818048925, 0.9118100898948334,
                      1.0110147142055477, 1.0180126386295232, 1.3201853093496474,
                      1.492172669340363]
        x = Analysis(df=failures_a, bounds='lrb', bounds_type='2s', cl=0.99)
        x.mle()

        # Exactly two points sharing a beta value is the normal case
        # (the eta-direction scan's own upper/lower solution pair).
        # More than two, spanning a wide eta range, is the quantization
        # artifact this test guards against.
        eta_span = x.eta_lrb.max() - x.eta_lrb.min()
        _, counts = np.unique(x.beta_lrb, return_counts=True)
        self.assertLessEqual(counts.max(), 4)

        for beta_val in np.unique(x.beta_lrb):
            etas_here = x.eta_lrb[x.beta_lrb == beta_val]
            if len(etas_here) > 2:
                self.assertLess((etas_here.max() - etas_here.min()) / eta_span, 0.01)


if __name__ == '__main__':
    unittest.main()
