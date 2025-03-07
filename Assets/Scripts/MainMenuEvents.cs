using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.UIElements;
public class MainMenuEvents : MonoBehaviour
{
    private UIDocument _document;
    private Button _button;
	public KeyCode Exit = KeyCode.Escape;

	private void Start()
    {
        _document = GetComponent<UIDocument>();

        _button = _document.rootVisualElement.Q("EnterButton") as Button;
	    _button.RegisterCallback<ClickEvent>(OnEnterButtonClick);

		_button = _document.rootVisualElement.Q("QuitButton") as Button;
		_button.RegisterCallback<ClickEvent>(OnQuitButtonClick);
	}
	private void Update()
	{
		if (Input.GetKeyDown(Exit))
		{
			if (_document.rootVisualElement.style.display == DisplayStyle.None)
			{
				_document.rootVisualElement.style.display = DisplayStyle.Flex;
			}
			else
			{
				_document.rootVisualElement.style.display = DisplayStyle.None;
			}
		}
	}

	private void OnEnterButtonClick(ClickEvent evt)
    {
		// Hides menu
		_document.rootVisualElement.style.display = DisplayStyle.None;
	}

	private void OnQuitButtonClick(ClickEvent evt)
	{
		// Closes application
		Application.Quit();
		Debug.Log("Application.Quit() called");

		
		//Temporary stop play mode in Unity editor
		#if UNITY_EDITOR
		UnityEditor.EditorApplication.isPlaying = false;
		#endif
	}
}