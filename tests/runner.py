import unittest

import tests.allow_clipping.test as allow_clipping_test
import tests.alternative_3_splice_site_novel.test as alternative_3_splice_site_novel_test
import tests.alternative_5_splice_site_novel.test as alternative_5_splice_site_novel_test
import tests.darts_model.test as darts_model_test
import tests.fixed_event_set.test as fixed_event_set_test
import tests.individual_counts.test as individual_counts_test
import tests.mutually_exclusive_exons_novel.test as mutually_exclusive_exons_novel_test
import tests.only_one_sample.test as only_one_sample_test
import tests.overlapped_gene.test as overlapped_gene_test
import tests.paired_stats.test as paired_stats_tests
import tests.prep_post.test as prep_post_test
import tests.read_count_edge_cases.test as read_count_edge_cases_test
import tests.retained_intron_novel.test as retained_intron_novel_test
import tests.skipped_exon_basic.test as skipped_exon_basic_test
import tests.skipped_exon_novel.test as skipped_exon_novel_test
import tests.stat_large_file.test as stat_large_file_test
import tests.stranded.test as stranded_test
import tests.task_stat.test as task_stat_test
import tests.variable_read_length.test as variable_read_length_test


def build_test_suite():
    loader = unittest.defaultTestLoader
    suite = unittest.TestSuite()
    suite.addTest(loader.loadTestsFromModule(allow_clipping_test))
    suite.addTest(
        loader.loadTestsFromModule(alternative_3_splice_site_novel_test))
    suite.addTest(
        loader.loadTestsFromModule(alternative_5_splice_site_novel_test))
    suite.addTest(loader.loadTestsFromModule(darts_model_test))
    suite.addTest(loader.loadTestsFromModule(fixed_event_set_test))
    suite.addTest(loader.loadTestsFromModule(individual_counts_test))
    suite.addTest(
        loader.loadTestsFromModule(mutually_exclusive_exons_novel_test))
    suite.addTest(loader.loadTestsFromModule(only_one_sample_test))
    suite.addTest(loader.loadTestsFromModule(overlapped_gene_test))
    suite.addTest(loader.loadTestsFromModule(paired_stats_tests))
    suite.addTest(loader.loadTestsFromModule(prep_post_test))
    suite.addTest(loader.loadTestsFromModule(read_count_edge_cases_test))
    suite.addTest(loader.loadTestsFromModule(retained_intron_novel_test))
    suite.addTest(loader.loadTestsFromModule(skipped_exon_basic_test))
    suite.addTest(loader.loadTestsFromModule(skipped_exon_novel_test))
    suite.addTest(loader.loadTestsFromModule(stat_large_file_test))
    suite.addTest(loader.loadTestsFromModule(stranded_test))
    suite.addTest(loader.loadTestsFromModule(task_stat_test))
    suite.addTest(loader.loadTestsFromModule(variable_read_length_test))
    return suite


def main():
    runner = unittest.TextTestRunner(verbosity=2)
    suite = build_test_suite()
    runner.run(suite)


if __name__ == '__main__':
    main()
