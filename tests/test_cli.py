from tacs_artifact_removal.cli import build_parser


def test_cli_keeps_time_domain_background_matching_opt_in():
    args = build_parser().parse_args(["recording.edf", "--frequency", "40"])

    assert args.background_match is False
    assert args.plot is False
    assert args.window == 1.0
    assert args.hop == 0.1
