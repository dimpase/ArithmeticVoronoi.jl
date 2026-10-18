"""Write readable presentations and a concise report from completed exact results."""
from fractions import Fraction
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def reduce_word(word):
    result = []
    for letter in word:
        if result and result[-1] == -letter:
            result.pop()
        else:
            result.append(letter)
    return result


def substitute(word, images):
    result = []
    for letter in word:
        image = images[abs(letter) - 1]
        result.extend(image if letter > 0 else [-x for x in reversed(image)])
    return reduce_word(result)


def word_text(word):
    terms = []
    i = 0
    while i < len(word):
        j = i + 1
        while j < len(word) and word[j] == word[i]:
            j += 1
        exponent = (j - i) * (1 if word[i] > 0 else -1)
        terms.append(f'g{abs(word[i])}' + ('' if exponent == 1 else f'^{exponent}'))
        i = j
    return ' '.join(terms) or '1'


def element_text(coefficients):
    a, b = map(Fraction, coefficients)
    if not b:
        return str(a)
    term = 'w' if abs(b) == 1 else f'{abs(b)}*w'
    first = str(a) if a else ''
    return first + (' - ' if b < 0 and first else '-' if b < 0 else ' + ' if first else '') + term


def main():
    lines = ['Arithmetic Voronoi: fields, constructive words and reusable data',
             '=' * 63, '', 'Completed arithmetic groups are K-linear GL2(O), SL2(O), and',
             'the separately constructed PSL2(O). GL/SL retain -I. The lattice is O^2.', '',
             'D   class h   perfect forms   oriented loops   GL small   SL small   PSL small']
    for D in (5, 8, 12, 13, 40):
        folder = ROOT / f'results/fields/D{D}'
        read = lambda name: json.loads((folder / name).read_text())
        run = read('run.json')
        assert all(run['status'][s] == 'complete' for s in ['enumeration', 'two_cells', 'GL', 'SL', 'PSL'])
        field = read('field.json')
        counts = run['presentations']
        size = lambda g: f'{counts[g]["simplified_generators"]}/{counts[g]["simplified_relators"]}'
        lines.append(f'{D:<3} {field["ordinary_class_number"]:<9} {run["geometry"]["perfect_form_orbits"]:<15}'
                     f'{run["geometry"]["oriented_loops"]:<17}{size("GL"):<11}{size("SL"):<11}{size("PSL")}')
        text = [f'Exact arithmetic group presentations, D={D}',
                f'w^2 - ({field["polynomial_coefficients"][1] * -1})*w + ({field["polynomial_coefficients"][0]}) = 0',
                'Matrices use (1,w); words are evaluated in written matrix-product order.',
                'GL and SL retain -I. PSL matrices are SL lifts modulo +/-I.', '']
        for name in ('GL', 'SL', 'PSL'):
            d = read(name + '_simplified.json')
            text.append(f'{name}: {len(d["generators"])} generators, {len(d["relators"])} relations')
            for g in d['generators']:
                rows = ['[' + ', '.join(element_text(x) for x in row) + ']' for row in g['matrix']]
                text.append(f'  g{g["id"]} = [' + '; '.join(rows) + ']')
            for i, w in enumerate(d['relators'], 1):
                text.append(f'  r{i}: {word_text(w)} = 1' + (' in PSL (lift may be +/-I)' if name == 'PSL' else ''))
            if name == 'SL':
                minus = substitute(read('SL_raw.json')['minus_identity_word'], d['old_generators_as_words_in_new'])
                text.append('  -I = ' + word_text(minus))
            text.append('')
        text += ['Both Tietze generator maps are in the corresponding *_simplified.json.',
                 'Completeness: ../../../PRESENTATION_JUSTIFICATION.md.',
                 'Stage statuses and timings: run.json; coverage and loops: two_cells/.',
                 'Constructive words: words/ and word_validation.json.',
                 'Kirchheimer-Wolfart comparison: KW/ (D5,D8,D12,D13 only).', '']
        (folder / 'PRESENTATIONS.txt').write_text('\n'.join(text))
    lines += ['', 'Small presentation sizes are generators / relators. Global two-cell orbit',
              'minimization is not performed: exhaustive local ridges in both orientations',
              'supply the complete redundant attaching loops.', '',
              'SL abelianizations (derived from the actual relators):']
    for D in (5, 8, 12, 13, 40):
        folder = ROOT / f'results/fields/D{D}'
        run = json.loads((folder / 'run.json').read_text())
        a = run['presentations']['SL']['abelianization']
        lines.append(f'  D{D}: free rank {a["free_rank"]}, torsion {a["torsion"]}')
    lines += ['', 'Recorded complete in-process pipeline timings (seconds):']
    for D in (5, 8, 12, 13, 40):
        folder = ROOT / f'results/fields/D{D}'
        initial = json.loads((folder / 'initial_complete_run.json').read_text())
        current = json.loads((folder / 'run.json').read_text())
        lines.append(f'  D{D}: initial complete run {initial["total_in_process_seconds"]:.3f}; '
                     f'enumeration cache={initial["enumeration_from_cache"]}; '
                     f'latest verified replay {current["total_in_process_seconds"]:.3f}.')
    lines += ['', 'D5/D8 initially reused preserved geometry. D12/D13 initial runs computed',
              'uncached enumeration and loops; first-use compilation and the original raw',
              'Smith-normal-form cost are included. D40 initial cold stages are separately',
              'retained in initial_cold_stage_checkpoint.json: 518.494 seconds for loops;',
              'raw and simplified presentations completed before interrupting an expensive',
              'raw abelianization. Its recorded complete run uses validated caches and',
              'the certified simplified presentation for the same invariant. There is no',
              'claimed uninterrupted cold end-to-end total for D40. All timings exclude',
              'Julia/OSCAR startup and the separate regression suite; one Julia thread.', '',
              'Words: exact general cone descent for valid matrices in these five fields,',
              'with GL/SL/PSL alphabets, resource errors, partial checkpoints, and exact',
              'complete-path cache replay. No canonical or shortest-word guarantee.', '',
              'Literature: explicit maps both ways for all four reference fields, proved',
              'by independently faithful presentations and commuting matrix maps. D8',
              'English (10.1) has -2 where the German transcription and exact relation',
              'give -5. No original-paper error is asserted without its scan.', '',
              'D40: ordinary class number two; the nonprincipal ideal (2,w) has norm two.',
              'The integer-primitive, O-nonunimodular vector (2,w) is included at trace 28',
              'by bounded lattice enumeration. No PID assumption enters the Voronoi engine.', '',
              'Verification: results/fields_delivery_validation.log, fields_words_tests.json,',
              'sanitizer_fields/summary.json and fields_dataset_complete_verification.json.',
              'Data: results/dataset/v1/dataset.json, with exact coefficients, all maps,',
              'certificates and the pinned checker under software/. LMFDB labels remain',
              'null pending verification; no submission, publication or contact was made.', '']
    (ROOT / 'results/FIELDS_WORDS_DATA_REPORT.txt').write_text('\n'.join(lines))
    print('\n'.join(lines[:13]))


if __name__ == '__main__':
    main()
