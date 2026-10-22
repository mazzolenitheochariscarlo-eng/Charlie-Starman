"""Genetic algorithm that breeds desk screen designs against a brief."""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .design import GENES, decode, evaluate


@dataclass
class Result:
    ranked: list  # Evaluations, best first, one per distinct design
    history: list  # (generation, best fitness, mean fitness)

    @property
    def best(self):
        return self.ranked[0]


def _tournament(rng, fitness, k=3):
    idx = rng.integers(0, len(fitness), size=k)
    return idx[np.argmax(fitness[idx])]


def run(brief, population=40, generations=30, seed=1, mutation_rate=0.15, elite=2, progress=None):
    rng = np.random.default_rng(seed)
    n_genes = len(GENES)
    pop = rng.uniform(size=(population, n_genes))
    archive = {}
    history = []

    def assess(vector):
        screen = decode(vector, brief)
        if screen not in archive:
            archive[screen] = evaluate(screen, brief)
        return archive[screen].fitness

    for gen in range(generations + 1):
        fitness = np.array([assess(v) for v in pop])
        history.append((gen, float(fitness.max()), float(fitness.mean())))
        if progress:
            progress(gen, fitness.max(), fitness.mean())
        if gen == generations:
            break
        order = np.argsort(fitness)[::-1]
        children = [pop[i].copy() for i in order[:elite]]
        while len(children) < population:
            a, b = pop[_tournament(rng, fitness)], pop[_tournament(rng, fitness)]
            mask = rng.uniform(size=n_genes) < 0.5
            child = np.where(mask, a, b)
            mutate = rng.uniform(size=n_genes) < mutation_rate
            child = np.where(mutate, np.clip(child + rng.normal(0, 0.15, n_genes), 0, 1), child)
            reset = rng.uniform(size=n_genes) < mutation_rate / 5  # occasional big jump
            child = np.where(reset, rng.uniform(size=n_genes), child)
            children.append(child)
        pop = np.array(children)

    ranked = sorted(archive.values(), key=lambda e: e.fitness, reverse=True)
    return Result(ranked, history)


def distinct(ranked, n=3):
    """Top designs that differ in pattern, felt or core, for a meaningful shortlist.

    Flat (non-formable) felts count as one pattern, since the cells then do nothing.
    """
    picked, seen = [], set()
    for e in ranked:
        pattern = e.screen.pattern if e.metrics["felt_formed"] else "flat"
        key = (pattern, e.screen.felt, e.screen.core)
        if key not in seen:
            seen.add(key)
            picked.append(e)
        if len(picked) == n:
            break
    return picked
