import testkd
import mock

def test_config():
    f = "config day 1 2"
    testkd.config(f)
    assert "day" == testkd.QUEUE
    assert "1" == testkd.MEM
    assert "2" == testkd.PROCESS


@mock.patch("testkd.pause")
@mock.patch("testkd.config")
@mock.patch("testkd.run_command")
def test_process_line_pause(mock_run_command, mock_config, mock_pause):
    testkd.process_line("pause")
    assert mock_pause.called
    assert not mock_run_command.called
    assert not mock_config.called


@mock.patch("testkd.pause")
@mock.patch("testkd.config")
@mock.patch("testkd.run_command")
def test_process_line_config(mock_run_command, mock_config, mock_pause):
    testkd.process_line("config 1 2 3")
    assert not mock_pause.called
    assert not mock_run_command.called
    assert mock_config.called


@mock.patch("testkd.pause")
@mock.patch("testkd.config")
@mock.patch("testkd.run_command")
def test_process_line_run(mock_run_command, mock_config, mock_pause):
    testkd.process_line("othercomand")
    assert not mock_pause.called
    assert mock_run_command.called
    assert not mock_config.called

